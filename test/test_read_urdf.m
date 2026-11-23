function test_read_urdf()
%TEST_READ_URDF Compare the XML reader with the Robotics System Toolbox reader.
repo = fileparts(fileparts(mfilename('fullpath')));
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
addpath(fullfile(repo, 'kinematics'), fullfile(repo, 'transformations'), fullfile(repo, 'dynamics'));
assert(~isempty(which('importrobot')), 'The comparison requires Robotics System Toolbox.');

urdf = fullfile(repo, 'urdf');
ur5e = fullfile(urdf, 'ur_description', 'urdf', 'ur5e.urdf');
dualArm = fullfile(urdf, 'lbr_description', 'urdf', 'iiwa7-dual-arm.urdf');
fixture = mixedJointURDF();
cases = {
    {ur5e}
    {string(ur5e), "tool0", "world"}
    {ur5e, 'tool0', 'shoulder_link'}
    {fullfile(urdf, 'ur_description', 'urdf', 'ur3e.urdf')}
    {fullfile(urdf, 'lbr_description', 'urdf', 'iiwa7.urdf')}
    {dualArm}
    {dualArm, 'left_link_ee', 'left_link_0'}
    {dualArm, 'right_link_ee', 'right_link_0'}
    {fullfile(urdf, 'IP', 'urdf', 'IP.urdf'), 'rod', 'base_link'}
    {fullfile(urdf, 'mass', 'urdf', 'mass.urdf')}
    {fullfile(urdf, 'two_link', 'urdf', 'two.urdf')}
    {fixture}
    {fixture, 'tool', 'payload'}
    {fixture, 'side'}
    };

expected = cell(size(cases));
for i = 1:numel(cases)
    args = cases{i};
    expected{i} = read_urdf_old(args{:});
    % Compare dynamics using the new reader's documented gravity default.
    expected{i}.gravity = [0, 0, -9.81];
    robot = read_urdf(args{:});
    assertRobotEqual(robot, expected{i}, i);
    q = linspace(-0.2, 0.3, robot.dof).';
    assertClose(forward_kin_general(robot, q), forward_kin_general(expected{i}, q), i, 'FK');
    assertClose(mass_matrix(robot, q), mass_matrix(expected{i}, q), i, 'mass matrix');
    if robot.dof > 0
        velocity = 0.4 * q;
        acceleration = -0.2 * q;
        wrench = (1:6).' / 10;
        assertClose(inverse_dynamics(robot, q, velocity, acceleration, wrench), ...
            inverse_dynamics(expected{i}, q, velocity, acceleration, wrench), i, 'inverse dynamics');
    end
end

invalid = {{fixture, 'missing'}, {fixture, 'tool', 'missing'}, ...
    {fixture, 'root'}, {fixture, 'tool', 'tool'}};
for i = 1:numel(invalid)
    previous = read_urdf_old(invalid{i}{:});
    assert(isempty(previous) && isequal(read_urdf(invalid{i}{:}), previous));
end
assertError(@() read_urdf_old(fixture, 'tool', 'side'));
assertError(@() read_urdf(fixture, 'tool', 'side'), 'read_urdf:InvalidChain');

% Run again with the toolbox off the MATLAB path, using the saved results.
entries = strsplit(path, pathsep);
roboticsDirectory = [fullfile(matlabroot, 'toolbox', 'robotics'), filesep];
roboticsEntries = entries(startsWith(entries, roboticsDirectory, 'IgnoreCase', ispc));
rmpath(roboticsEntries{:});
assert(isempty(which('importrobot')) && isempty(which('rigidBodyTree')));
for i = 1:numel(cases)
    assertRobotEqual(read_urdf(cases{i}{:}), expected{i}, i);
end
fprintf('read_urdf: %d valid cases and %d invalid endpoint cases matched the legacy reader.\n', ...
    numel(cases), numel(invalid) + 1);
end

function assertError(call, identifier)
try
    call();
catch err
    if nargin >= 2
        assert(strcmp(err.identifier, identifier), ...
            'Expected error %s, received %s.', identifier, err.identifier);
    end
    return;
end
error('test_read_urdf:MissingError', 'Expected the reader to reject an invalid chain.');
end

function assertRobotEqual(actual, expected, index)
assert(isstruct(actual) && isstruct(expected) && isequal(fieldnames(actual), fieldnames(expected)), ...
    'Case %d: unexpected robot fields.', index);
fields = fieldnames(expected);
for i = 1:numel(fields)
    name = fields{i};
    assertClose(actual.(name), expected.(name), index, name);
end
end

function assertClose(actual, expected, index, field)
tolerance = 1e-12 + 1e-12 * abs(expected(:));
assert(isequal(size(actual), size(expected)) && all(abs(actual(:) - expected(:)) <= tolerance), ...
    'Case %d: %s differs from the legacy reader.', index, field);
end

function xml = mixedJointURDF()
xml = strjoin({
    '<robot name="mixed_joint_test">'
    '  <vendor:metadata/>'
    '  <link name="tool"/>'
    '  <link name="root"/>'
    '  <link name="moving"/>'
    '  <link name="payload">'
    '    <inertial>'
    '      <origin xyz="0.1 -0.2 0.3" rpy="0.2 0.3 -0.4"/>'
    '      <mass value="2.1"/>'
    '      <inertia ixx="0.5" iyy="0.6" izz="0.7" ixy="0.02" ixz="0.03" iyz="-0.04"/>'
    '    </inertial>'
    '  </link>'
    '  <link name="slide"/>'
    '  <link name="side"/>'
    '  <joint name="z_branch" type="continuous">'
    '    <parent link="root"/><child link="moving"/>'
    '    <origin xyz="0 0 0.5" rpy="0 0.2 0"/>'
    '  </joint>'
    '  <joint name="mount" type="fixed">'
    '    <parent link="moving"/><child link="payload"/>'
    '    <origin xyz="0 0.4 0" rpy="0.1 0.2 0.3"/>'
    '  </joint>'
    '  <joint name="slider" type="prismatic">'
    '    <parent link="payload"/><child link="slide"/>'
    '    <axis xyz="0 -3 4"/><limit lower="-1" upper="1" effort="10" velocity="2"/>'
    '  </joint>'
    '  <joint name="tcp" type="fixed">'
    '    <parent link="slide"/><child link="tool"/><origin xyz="0.1 0 0.2"/>'
    '  </joint>'
    '  <joint name="a_branch" type="fixed">'
    '    <parent link="root"/><child link="side"/><origin xyz="0 0 -1"/>'
    '  </joint>'
    '</robot>'
    }, newline);
end
