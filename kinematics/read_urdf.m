function robot = read_urdf(filename, flange_link, base_link)
%READ_URDF Build a robot chain from a URDF file or URDF XML text.
%   The base defaults to the root link. The flange defaults to the last
%   link in depth-first order, with siblings ordered by joint name.
%   The chain supports fixed, revolute, continuous, and prismatic joints.
%   Gravity in the base frame is [0, 0, -9.81] m/s^2.

robot = [];
validateattributes(filename, {'char', 'string'}, {'scalartext', 'nonempty'});
source = char(filename);
parser = matlab.io.xml.dom.Parser;
% Older URDFs may have extension tags with undeclared namespace prefixes.
parser.Configuration.Namespaces = false;
if startsWith(strtrim(source), '<')
    document = parseString(parser, source);
else
    if ~isfile(source)
        source = which(source);
    end
    if isempty(source)
        error('read_urdf:FileNotFound', 'URDF file not found: %s', filename);
    end
    [~, ~, extension] = fileparts(source);
    if strcmpi(extension, '.xacro')
        error('read_urdf:UnexpandedXacro', 'Expand the Xacro file to URDF before reading it.');
    end
    document = parseFile(parser, source);
end
root = document.getDocumentElement();
if ~strcmp(root.getTagName(), 'robot') || isempty(root.getAttribute('name'))
    error('read_urdf:InvalidURDF', 'Expected a named URDF <robot> root element.');
end

linkNodes = childElements(root, 'link');
jointNodes = childElements(root, 'joint');
linkNames = elementNames(linkNodes, 'link');
jointNames = elementNames(jointNodes, 'joint');
[parents, incomingJoint, bodyOrder] = readTree(linkNames, jointNodes, jointNames);
rootLink = bodyOrder(1);

base = rootLink;
if nargin >= 3
    base = find(strcmp(linkNames, base_link), 1);
    if isempty(base)
        disp('not a valid base');
        return;
    end
end
flange = bodyOrder(end);
if nargin >= 2
    flange = find(strcmp(linkNames, flange_link), 1);
    if isempty(flange) || flange == rootLink
        disp('not a valid flange');
        return;
    end
elseif flange == rootLink
    error('read_urdf:InvalidChain', 'The robot has no link that can be used as the flange.');
end

chain = zeros(1, numel(linkNodes));
head = numel(chain) + 1;
link = flange;
while link ~= base && link ~= rootLink
    head = head - 1;
    chain(head) = link;
    link = parents(link);
end
chain = chain(head:end);
if link ~= base
    error('read_urdf:InvalidChain', 'Flange "%s" is not a descendant of base "%s".', ...
        linkNames{flange}, linkNames{base});
end
if isempty(chain)
    disp('not a valid chain');
    return;
end

bodyList = cell(1, numel(chain));
for i = 1:numel(chain)
    link = chain(i);
    bodyList{i} = readBody(linkNodes{link}, jointNodes{incomingJoint(link)});
end
robot = build_robot_chain(bodyList);
robot.gravity = [0, 0, -9.81];
end

function [parents, incomingJoint, order] = readTree(linkNames, joints, jointNames)
n = numel(linkNames);
parents = zeros(1, n);
incomingJoint = zeros(1, n);
children = cell(1, n);
[~, jointOrder] = sort(jointNames);
for i = jointOrder
    joint = joints{i};
    parent = find(strcmp(linkNames, attribute(firstChild(joint, 'parent'), 'link')), 1);
    child = find(strcmp(linkNames, attribute(firstChild(joint, 'child'), 'link')), 1);
    if isempty(parent) || isempty(child) || parent == child || parents(child) ~= 0
        error('read_urdf:InvalidTree', 'Invalid parent or child link for joint "%s".', jointNames{i});
    end
    parents(child) = parent;
    incomingJoint(child) = i;
    children{parent}(end + 1) = child;
end

root = find(parents == 0);
if numel(root) ~= 1
    error('read_urdf:InvalidTree', 'A URDF robot must have one root link.');
end
order = zeros(1, n);
count = 0;
stack = root;
while ~isempty(stack)
    link = stack(end);
    stack(end) = [];
    count = count + 1;
    order(count) = link;
    stack = [stack, fliplr(children{link})]; %#ok<AGROW>
end
if count ~= n
    error('read_urdf:InvalidTree', 'Every link must be connected to the root link.');
end
end

function body = readBody(link, jointNode)
joint = struct('Name', jointNode.getAttribute('name'), ...
    'Type', jointNode.getAttribute('type'), 'JointAxis', zeros(1, 3), ...
    'JointToParentTransform', originTransform(jointNode), 'ChildToJointTransform', eye(4));
switch joint.Type
    case {'revolute', 'continuous', 'prismatic'}
        axis = numericAttribute(firstChild(jointNode, 'axis'), 'xyz', [1 0 0]);
        lengthAxis = norm(axis);
        if lengthAxis == 0 || ~isfinite(lengthAxis)
            error('read_urdf:InvalidAxis', 'Invalid axis for joint "%s".', joint.Name);
        end
        joint.JointAxis = axis / lengthAxis;
        if strcmp(joint.Type, 'continuous')
            joint.Type = 'revolute';
        end
    case 'fixed'
        % A fixed joint contributes only a transform and inertial properties.
    otherwise
        error('read_urdf:UnsupportedJoint', 'Unsupported joint type "%s" for "%s".', ...
            joint.Type, joint.Name);
end

body = struct('Joint', joint, 'Mass', 0, 'CenterOfMass', zeros(1, 3), 'Inertia', zeros(1, 6));
inertial = firstChild(link, 'inertial');
if isempty(inertial)
    return;
end
body.Mass = numericAttribute(firstChild(inertial, 'mass'), 'value');
if body.Mass < 0
    error('read_urdf:InvalidMass', 'Negative mass for link "%s".', link.getAttribute('name'));
end
pose = originTransform(inertial);
body.CenterOfMass = pose(1:3, 4).';

inertia = firstChild(inertial, 'inertia');
ixx = numericAttribute(inertia, 'ixx');
iyy = numericAttribute(inertia, 'iyy');
izz = numericAttribute(inertia, 'izz');
ixy = numericAttribute(inertia, 'ixy');
ixz = numericAttribute(inertia, 'ixz');
iyz = numericAttribute(inertia, 'iyz');
Icom = [ixx ixy ixz; ixy iyy iyz; ixz iyz izz];
% URDF specifies inertia at the center of mass. The chain uses inertia
% about the link origin, expressed in the link frame.
I = transform_inertia_matrix(Icom, body.Mass, pose);
body.Inertia = [I(1,1), I(2,2), I(3,3), I(2,3), I(1,3), I(1,2)];
end

function pose = originTransform(node)
origin = firstChild(node, 'origin');
xyz = numericAttribute(origin, 'xyz', [0 0 0]);
rpy = numericAttribute(origin, 'rpy', [0 0 0]);
pose = make_tform(RPY2R(rpy), xyz);
end

function names = elementNames(elements, tag)
names = cellfun(@(element) element.getAttribute('name'), elements, 'UniformOutput', false);
if any(cellfun(@isempty, names)) || numel(unique(names)) ~= numel(names)
    error('read_urdf:InvalidName', 'Each <%s> must have a unique name.', tag);
end
end

function elements = childElements(node, tag)
elements = {};
child = node.getFirstElementChild();
while ~isempty(child)
    if strcmp(child.getTagName(), tag)
        elements{end + 1} = child; %#ok<AGROW>
    end
    child = child.getNextElementSibling();
end
end

function child = firstChild(node, tag)
child = node.getFirstElementChild();
while ~isempty(child) && ~strcmp(child.getTagName(), tag)
    child = child.getNextElementSibling();
end
end

function value = attribute(node, name)
value = '';
if ~isempty(node)
    value = node.getAttribute(name);
end
end

function values = numericAttribute(node, name, default)
text = strtrim(attribute(node, name));
if nargin < 3
    count = 1;
else
    count = numel(default);
    if isempty(text)
        values = default;
        return;
    end
end
values = str2double(regexp(text, '\s+', 'split'));
if numel(values) ~= count || ~isreal(values) || any(~isfinite(values))
    error('read_urdf:InvalidAttribute', 'Expected %d finite value(s) for "%s".', count, name);
end
end
