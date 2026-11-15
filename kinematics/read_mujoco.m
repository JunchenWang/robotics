function [robot, joint_names] = read_mujoco(filename, flange_link, base_link)
%READ_MUJOCO Build a robot chain from an MJCF XML file using MATLAB.
%   Body or site names select the flange and base. Defaults are the last
%   body and 'world'. If a body and site share a name, use 'site:NAME'
%   to select the site. Joint coordinates are offsets from MJCF ref.
%   Supports hinge and slide joints, explicit inertials, and inertia from
%   solid sphere, ellipsoid, box, cylinder, and capsule geoms.
%   Geoms with mesh or shell inertia require an explicit <inertial>.

robot = [];
joint_names = {};
if ~isfile(filename)
    error('read_mujoco:FileNotFound', 'XML file not found: %s', filename);
end
[~, fileInfo] = fileattrib(filename);
parser = matlab.io.xml.dom.Parser;
doc = parseFile(parser, fileInfo.Name);
documentRoot = doc.getDocumentElement();
if ~strcmp(documentRoot.getTagName(), 'mujoco')
    error('read_mujoco:InvalidXML', 'Expected an MJCF <mujoco> root element.');
end
root = readElement(documentRoot, parser, fileparts(fileInfo.Name), {fileInfo.Name});

compiler = struct();
for node = elements(root, 'compiler')
    compiler = attributes(node{1}, compiler);
end
if strcmp(attribute(compiler, 'coordinate', 'local'), 'global')
    error('read_mujoco:UnsupportedCompiler', 'Global MJCF coordinates are unsupported.');
end
for field = {'boundmass', 'boundinertia', 'settotalmass'}
    if numbers(compiler, field{1}, 0) > 0
        error('read_mujoco:UnsupportedCompiler', ...
            'Compiler attribute "%s" changes model inertias and is unsupported.', field{1});
    end
end
if strcmp(attribute(compiler, 'balanceinertia', 'false'), 'true')
    error('read_mujoco:UnsupportedCompiler', 'Compiler balanceinertia is unsupported.');
end

defaults = containers.Map('KeyType', 'char', 'ValueType', 'any');
defaults('main') = struct('joint', struct(), 'geom', struct(), 'site', struct());
for node = elements(root, 'default')
    loadDefaults(node{1}, defaults, defaults('main'), 'main');
end

bodies = struct('name', 'world', 'parent', 0, 'pose', eye(4), ...
    'class', 'main', 'node', []);
sites = struct('name', {}, 'parent', {}, 'pose', {}, 'class', {}, 'props', {});
for node = elements(root, 'worldbody')
    world = node{1};
    [bodies, sites] = readBodies(world, bodies, sites, 1, eye(4), ...
        attribute(world.props, 'childclass', 'main'), compiler);
end
base = 1;
basePose = eye(4);
baseIsSite = false;
if nargin >= 3
    [base, basePose, baseIsSite] = resolveEndpoint(base_link, bodies, sites, compiler, defaults);
    if isempty(base)
        warning('read_mujoco:InvalidBase', 'Unknown base body or site: %s', base_link);
        return;
    end
end
flange = numel(bodies);
flangePose = eye(4);
flangeIsSite = false;
if nargin >= 2
    [flange, flangePose, flangeIsSite] = resolveEndpoint( ...
        flange_link, bodies, sites, compiler, defaults);
    if isempty(flange)
        warning('read_mujoco:InvalidFlange', 'Unknown flange body or site: %s', flange_link);
        return;
    end
end

chain = zeros(1, numel(bodies));
head = numel(chain) + 1;
body = flange;
while body ~= base && body ~= 1
    head = head - 1;
    chain(head) = body;
    body = bodies(body).parent;
end
chain = chain(head:end);
if body ~= base || (isempty(chain) && ~baseIsSite && ~flangeIsSite)
    warning('read_mujoco:InvalidChain', ...
        'The flange must be attached to the base body or one of its descendants.');
    return;
end

% Express the world gravity in the selected base frame at the reference pose.
gravity = [0 0 -9.81];
gravityEnabled = true;
for node = elements(root, 'option')
    gravity = numbers(attributes(node{1}), 'gravity', gravity);
    for flag = elements(node{1}, 'flag')
        state = attribute(flag{1}.props, 'gravity', '');
        if ~isempty(state)
            gravityEnabled = strcmp(state, 'enable');
        end
    end
end
worldRotation = basePose(1:3, 1:3);
body = base;
while body ~= 1
    worldRotation = bodies(body).pose(1:3, 1:3) * worldRotation;
    body = bodies(body).parent;
end

robot.dof = 0;
robot.mass = zeros(0, 1);
robot.inertia = zeros(3, 3, 0);
robot.A = zeros(0, 6);
robot.M = zeros(4, 4, 0);
robot.com = zeros(0, 3);
robot.gravity = gravityEnabled * gravity * worldRotation;
robot.TCP = eye(4);

link = 0;
offset = tform_inv(basePose);
for body = bodies(chain)
    offset = offset * body.pose;
    parts = bodyParts(body.node, eye(4), body.class, compiler);
    for part = parts(strcmp({parts.tag}, 'joint') | strcmp({parts.tag}, 'freejoint'))
        props = withDefaults(part, defaults, 'joint');
        type = attribute(props, 'type', 'hinge');
        if strcmp(part.tag, 'freejoint')
            type = 'free';
        end
        axis = part.pose(1:3, 1:3) * unitVector(numbers(props, 'axis', [0 0 1])).';
        switch type
            case 'hinge'
                anchor = part.pose * [numbers(props, 'pos', [0 0 0]).'; 1];
                screw = [axis; cross(anchor(1:3), axis)];
            case 'slide'
                screw = [zeros(3, 1); axis];
            otherwise
                error('read_mujoco:UnsupportedJoint', ...
                    'Joint type "%s" on body "%s" is unsupported.', type, body.name);
        end
        % All joints on a body share its reference frame. The body's mass
        % and inertia belong to the last joint because it moves with all of them.
        link = link + 1;
        robot.mass(link, 1) = 0;
        robot.inertia(:, :, link) = zeros(3);
        robot.com(link, :) = zeros(1, 3);
        joint_names{link} = attribute(props, 'name', '');
        robot.A(link, :) = screw.';
        robot.M(:, :, link) = offset;
        offset = eye(4);
    end
    if link > 0
        [mass, firstMoment, inertia] = bodyInertia(parts, offset, compiler, defaults);
        totalMass = robot.mass(link) + mass;
        if totalMass > 0
            robot.com(link, :) = (robot.mass(link) * robot.com(link, :) + ...
                firstMoment) / totalMass;
        end
        robot.mass(link) = totalMass;
        robot.inertia(:, :, link) = robot.inertia(:, :, link) + inertia;
    end
end
robot.dof = link;
robot.ME = offset * flangePose;
end

function [body, pose, isSite] = resolveEndpoint(name, bodies, sites, compiler, defaults)
name = char(name);
pose = eye(4);
body = [];
isSite = false;
if ~startsWith(name, 'site:')
    if isempty(name)
        return;
    end
    body = find(strcmp({bodies.name}, name), 1);
    if ~isempty(body)
        return;
    end
else
    name = name(6:end);
end
if isempty(name)
    return;
end
index = find(strcmp({sites.name}, name), 1);
if ~isempty(index)
    site = sites(index);
    body = site.parent;
    props = withDefaults(site, defaults, 'site');
    if isfield(props, 'fromto')
        [localPose, ~] = fromToPose(props);
    else
        localPose = readPose(props, compiler);
    end
    pose = site.pose * localPose;
    isSite = true;
end
end

function node = readElement(element, parser, modelDir, stack)
node = struct('tag', element.getTagName(), 'props', struct(), 'children', {{}});
items = element.getAttributes();
for i = 0:items.getLength() - 1
    item = items.item(i);
    name = item.getNodeName();
    if isvarname(name)
        node.props.(name) = item.getNodeValue();
    end
end
child = element.getFirstElementChild();
while ~isempty(child)
    if strcmp(child.getTagName(), 'include')
        filename = child.getAttribute('file');
        if ~isempty(regexp(filename, '^(?:[A-Za-z]:|[/\\])', 'once'))
            includePath = filename;
        else
            includePath = fullfile(modelDir, filename);
        end
        [exists, info] = fileattrib(includePath);
        if ~exists || any(strcmp(stack, info.Name))
            error('read_mujoco:InvalidInclude', 'Missing or recursive include: %s', filename);
        end
        included = parseFile(parser, info.Name);
        content = readElement(included.getDocumentElement(), parser, modelDir, [stack, {info.Name}]);
        node.children = [node.children, content.children];
    else
        node.children{end + 1} = readElement(child, parser, modelDir, stack);
    end
    child = child.getNextElementSibling();
end
end

function loadDefaults(node, defaults, inherited, fallbackName)
class = attribute(attributes(node), 'class', fallbackName);
for tag = {'joint', 'geom', 'site'}
    key = tag{1};
    for child = elements(node, key)
        inherited.(key) = attributes(child{1}, inherited.(key));
    end
end
defaults(class) = inherited;
for child = elements(node, 'default')
    loadDefaults(child{1}, defaults, inherited, '');
end
end

function [bodies, sites] = readBodies(node, bodies, sites, parent, frame, class, compiler)
for child = elements(node)
    child = child{1};
    props = attributes(child);
    childClass = attribute(props, 'childclass', class);
    switch child.tag
        case 'body'
            pose = frame * readPose(props, compiler);
            bodies(end + 1) = struct('name', attribute(props, 'name', ''), ...
                'parent', parent, 'pose', pose, 'class', childClass, 'node', child);
            [bodies, sites] = readBodies(child, bodies, sites, numel(bodies), ...
                eye(4), childClass, compiler);
        case 'frame'
            pose = frame * readPose(props, compiler);
            [bodies, sites] = readBodies(child, bodies, sites, parent, ...
                pose, childClass, compiler);
        case 'site'
            sites(end + 1) = struct('name', attribute(props, 'name', ''), ...
                'parent', parent, 'pose', frame, 'class', class, 'props', props);
        case {'attach', 'replicate', 'flexcomp', 'composite'}
            error('read_mujoco:UnsupportedElement', ...
                'Model generation element <%s> is unsupported.', child.tag);
    end
end
end

function parts = bodyParts(node, frame, class, compiler)
parts = struct('tag', {}, 'props', {}, 'pose', {}, 'class', {});
for child = elements(node)
    tag = child{1}.tag;
    props = attributes(child{1});
    switch tag
        case 'frame'
            parts = [parts, bodyParts(child{1}, frame * readPose(props, compiler), ...
                attribute(props, 'childclass', class), compiler)];
        case {'joint', 'freejoint', 'geom', 'inertial'}
            parts(end + 1) = struct('tag', tag, 'props', props, 'pose', frame, 'class', class);
    end
end
end

function props = withDefaults(part, defaults, tag)
class = attribute(part.props, 'class', part.class);
if ~isKey(defaults, class)
    error('read_mujoco:UnknownClass', 'Unknown MJCF defaults class: %s', class);
end
preset = defaults(class);
props = overlay(preset.(tag), part.props);
end

function [mass, firstMoment, inertia] = bodyInertia(parts, offset, compiler, defaults)
mass = 0;
firstMoment = zeros(1, 3);
inertia = zeros(3);
explicit = parts(strcmp({parts.tag}, 'inertial'));
mode = attribute(compiler, 'inertiafromgeom', 'auto');
if ~isempty(explicit) && ~strcmp(mode, 'true')
    props = explicit.props;
    m = numbers(props, 'mass', 0);
    pose = offset * explicit.pose * readPose(props, compiler);
    if isfield(props, 'fullinertia')
        v = numbers(props, 'fullinertia', zeros(1, 6));
        central = [v(1), v(4), v(5); v(4), v(2), v(6); v(5), v(6), v(3)];
    else
        central = diag(numbers(props, 'diaginertia', zeros(1, 3)));
    end
    mass = m;
    firstMoment = m * pose(1:3, 4).';
    inertia = transform_inertia_matrix(central, m, pose);
elseif ~strcmp(mode, 'false')
    groups = numbers(compiler, 'inertiagrouprange', [0 5]);
    for part = parts(strcmp({parts.tag}, 'geom'))
        props = withDefaults(part, defaults, 'geom');
        group = numbers(props, 'group', 0);
        if group >= groups(1) && group <= groups(2)
            [m, central, pose] = geomInertia(props, compiler);
            pose = offset * part.pose * pose;
            mass = mass + m;
            firstMoment = firstMoment + m * pose(1:3, 4).';
            inertia = inertia + transform_inertia_matrix(central, m, pose);
        end
    end
end
end

function [mass, inertia, pose] = geomInertia(props, compiler)
type = attribute(props, 'type', 'sphere');
sizeValues = numbers(props, 'size', [], 0);
sizeValues(end+1:3) = 0;
r = sizeValues(1);
h = sizeValues(2);
mass = numbers(props, 'mass', []);
density = numbers(props, 'density', 1000);
inertia = zeros(3);
pose = eye(4);
if isequal(mass, 0) || (isempty(mass) && density == 0) || ...
        any(strcmp(type, {'plane', 'hfield'}))
    mass = 0;
    return;
end
if isfield(props, 'mesh') || strcmp(attribute(props, 'shellinertia', 'false'), 'true')
    error('read_mujoco:UnsupportedInertia', ...
        'Specify an <inertial> for bodies that use mesh or shell inertia.');
end
if isfield(props, 'fromto')
    [pose, h] = fromToPose(props);
    if any(strcmp(type, {'ellipsoid', 'box'}))
        sizeValues = [r, r, h];
    end
else
    pose = readPose(props, compiler);
end
switch type
    case 'sphere'
        volume = 4*pi*r^3/3;
        moment = repmat(2*r^2/5, 1, 3);
    case {'box', 'ellipsoid'}
        moment = [sum(sizeValues(2:3).^2), r^2 + sizeValues(3)^2, ...
            sum(sizeValues(1:2).^2)];
        if strcmp(type, 'box')
            volume = 8*prod(sizeValues);
            moment = moment / 3;
        else
            volume = 4*pi*prod(sizeValues)/3;
            moment = moment / 5;
        end
    case 'cylinder'
        volume = 2*pi*r^2*h;
        moment = [r^2/4 + h^2/3, r^2/4 + h^2/3, r^2/2];
    case 'capsule'
        cylinderVolume = 2*pi*r^2*h;
        sphereVolume = 4*pi*r^3/3;
        volume = cylinderVolume + sphereVolume;
        axial = (cylinderVolume*r^2/2 + sphereVolume*2*r^2/5)/volume;
        transverse = (cylinderVolume*(r^2/4 + h^2/3) + ...
            sphereVolume*(2*r^2/5 + h^2 + 3*h*r/4))/volume;
        moment = [transverse, transverse, axial];
    otherwise
        error('read_mujoco:UnsupportedInertia', ...
            'Specify an <inertial> for bodies that use geom type "%s".', type);
end
if volume <= 0
    error('read_mujoco:InvalidGeom', 'Geom "%s" must have positive size.', type);
end
if isempty(mass)
    mass = density * volume;
end
inertia = diag(mass * moment);
end

function [pose, halfLength] = fromToPose(props)
endpoints = numbers(props, 'fromto', zeros(1, 6));
direction = endpoints(4:6) - endpoints(1:3);
center = (endpoints(4:6) + endpoints(1:3)) / 2;
halfLength = norm(direction) / 2;
pose = make_tform(zRotation(direction), center);
end

function pose = readPose(props, compiler)
angleScale = 1;
if strcmp(attribute(compiler, 'angle', 'degree'), 'degree')
    angleScale = pi / 180;
end
R = eye(3);
if isfield(props, 'quat')
    R = RMByQuaternion(unitVector(numbers(props, 'quat', [1 0 0 0])));
elseif isfield(props, 'axisangle')
    v = numbers(props, 'axisangle', zeros(1, 4));
    R = exp_w(unitVector(v(1:3)) * v(4) * angleScale);
elseif isfield(props, 'euler')
    angles = numbers(props, 'euler', zeros(1, 3)) * angleScale;
    sequence = attribute(compiler, 'eulerseq', 'xyz');
    directions = eye(3);
    for i = 1:3
        axis = directions(lower(sequence(i)) - 'x' + 1, :);
        rotation = exp_w(axis * angles(i));
        if isstrprop(sequence(i), 'lower')
            R = R * rotation;
        else
            R = rotation * R;
        end
    end
elseif isfield(props, 'xyaxes')
    v = numbers(props, 'xyaxes', zeros(1, 6));
    x = unitVector(v(1:3));
    z = unitVector(cross(x, v(4:6)));
    R = [x; cross(z, x); z].';
elseif isfield(props, 'zaxis')
    R = zRotation(numbers(props, 'zaxis', [0 0 1]));
end
pose = make_tform(R, numbers(props, 'pos', [0 0 0]));
end

function R = zRotation(direction)
z = unitVector(direction);
axis = cross([0 0 1], z);
if norm(axis) <= eps
    R = diag([1, sign(z(3)), sign(z(3))]);
else
    R = exp_w(unitVector(axis) * atan2(norm(axis), z(3)));
end
end

function v = unitVector(v)
length = norm(v);
if length <= eps
    error('read_mujoco:InvalidDirection', 'Expected a nonzero axis or quaternion.');
end
v = v / length;
end

function list = elements(node, tag)
list = node.children;
if nargin >= 2
    list = list(cellfun(@(child) strcmp(child.tag, tag), list));
end
end

function props = attributes(node, inherited)
if nargin < 2
    inherited = struct();
end
props = overlay(inherited, node.props);
end

function props = overlay(props, actual)
orientations = {'quat', 'axisangle', 'euler', 'xyaxes', 'zaxis'};
if any(isfield(actual, orientations))
    props = rmfield(props, intersect(fieldnames(props), orientations));
end
for key = fieldnames(actual).'
    props.(key{1}) = actual.(key{1});
end
end

function value = attribute(props, key, fallback)
value = fallback;
if isfield(props, key)
    value = props.(key);
end
end

function value = numbers(props, key, fallback, count)
if nargin < 4
    count = max(1, numel(fallback));
end
value = fallback;
if isfield(props, key)
    tokens = regexp(strtrim(props.(key)), '\s+', 'split');
    value = str2double(tokens);
    if any(~isfinite(value)) || (count > 0 && numel(value) ~= count)
        error('read_mujoco:InvalidAttribute', 'Invalid numeric attribute "%s".', key);
    end
end
end
