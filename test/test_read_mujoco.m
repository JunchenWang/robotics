function test_read_mujoco(outputDir)
%TEST_READ_MUJOCO Compare MJCF, exported URDF, and the MuJoCo runtime.
%   Python with the mujoco package must be available on PATH. Pass an
%   output directory to keep the converted URDFs and reference values.

testDir = fileparts(mfilename('fullpath'));
repo = fileparts(testDir);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
addpath(fullfile(repo, 'kinematics'), fullfile(repo, 'transformations'), ...
    fullfile(repo, 'quaternions'), fullfile(repo, 'dynamics'));
if nargin == 0
    outputDir = tempname;
    mkdir(outputDir);
    outputCleanup = onCleanup(@() rmdir(outputDir, 's')); %#ok<NASGU>
end
script = fullfile(testDir, 'read_mujoco_reference.py');
[status, output] = system(sprintf('python "%s" "%s"', script, outputDir));
assert(status == 0, '%s', output);
fprintf('%s', output);
cases = jsondecode(fileread(fullfile(outputDir, 'reference.json')));
for k = 1:numel(cases)
    sampleCase = cases(k);
    [mjcf, names] = read_mujoco(sampleCase.mjcf, sampleCase.flange, sampleCase.base);
    urdf = read_urdf(sampleCase.urdf, sampleCase.urdf_flange, sampleCase.urdf_base);
    label = sampleCase.label;
    assert(isequal(names(:), sampleCase.joint_names(:)), '%s: joint names or order differ.', label);
    assertClose(mjcf.gravity, sampleCase.gravity(:).', label, 'gravity');

    % The repository's URDF reader supplies a fixed gravity default. Use
    % the same base-frame gravity for both rigid-body dynamics calculations.
    urdf.gravity = sampleCase.gravity(:).';
    % URDF active link origins are at their joint anchors. Re-express the
    % MJCF struct in those frames before comparing the stored parameters.
    reframed = atJointAnchors(mjcf, reshape(sampleCase.anchors, 3, []).');
    fieldError = 0;
    for field = fieldnames(urdf).'
        name = field{1};
        fieldError = max(fieldError, assertClose(reframed.(name), urdf.(name), label, name));
    end

    errors = zeros(1, 3);
    for i = 1:numel(sampleCase.samples)
        sample = sampleCase.samples(i);
        q = sample.q(:);
        qd = sample.qd(:);
        qdd = sample.qdd(:);
        expectedPose = reshape(sample.pose, 4, 4);
        expectedMass = reshape(sample.mass, mjcf.dof, mjcf.dof);
        for imported = {mjcf, urdf}
            robot = imported{1};
            pose = forward_kin_general(robot, q);
            mass = mass_matrix(robot, q);
            force = inverse_dynamics(robot, q, qd, qdd, zeros(6, 1));
            errors = max(errors, [assertClose(pose, expectedPose, label, 'pose'), ...
                assertClose(mass, expectedMass, label, 'mass matrix'), ...
                assertClose(force, sample.force(:), label, 'inverse dynamics')]);
        end
    end
    fprintf('%s (%d DOF, %d states): fields %.3g, pose %.3g, mass matrix %.3g, inverse dynamics %.3g\n', ...
        label, mjcf.dof, numel(sampleCase.samples), fieldError, errors);
end
end

function robot = atJointAnchors(robot, anchors)
previous = eye(4);
for i = 1:robot.dof
    offset = make_tform(eye(3), anchors(i, :));
    robot.M(:, :, i) = tform_inv(previous) * robot.M(:, :, i) * offset;
    robot.A(i, :) = (adjoint_T(tform_inv(offset)) * robot.A(i, :).').';
    m = robot.mass(i);
    if m > 0
        inertia = transform_com_inertia_matrix(robot.inertia(:, :, i), m, ...
            make_tform(eye(3), robot.com(i, :)));
        robot.com(i, :) = robot.com(i, :) - anchors(i, :);
        robot.inertia(:, :, i) = transform_inertia_matrix(inertia, m, ...
            make_tform(eye(3), robot.com(i, :)));
    end
    previous = offset;
end
robot.ME = tform_inv(previous) * robot.ME;
end

function maxError = assertClose(actual, expected, label, field)
assert(isequal(size(actual), size(expected)), '%s: unexpected size of %s.', label, field);
maxError = max(abs(actual(:) - expected(:)));
% MuJoCo diagonalizes fullinertia and reconstructing it introduces roundoff.
bound = 1e-9 + 1e-10 * max(abs(expected(:)));
assert(maxError <= bound, '%s: %s maximum error %.16g exceeds %.3g.', label, field, maxError, bound);
end
