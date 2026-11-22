function pose_ik_feedback_step(block)
%POSE_IK_FEEDBACK_STEP Discrete Level-2 MATLAB S-function for task-space position IK.
% Dialog parameters: robot (from read_robot), cfg (see mj_ur5e_pose_ik_feedback.slx).
% Inputs: 4-by-4 desired tool pose in the base frame and measured joints (rad).
% Outputs: joint position command (rad) and measured-pose convergence.
% Each sample takes one damped least-squares step from the measured joints.
% cfg.qd_max optionally overrides the UR5e default of pi rad/s. Use a
% scalar or one value per joint (rad/s for revolute, m/s for prismatic).
% Use +Inf to disable speed limiting when the velocity controller's input
% is limited by the caller.

block.NumDialogPrms = 2;
block.DialogPrmsTunable = {'Nontunable', 'Nontunable'};
robot = block.DialogPrm(1).Data;
cfg = block.DialogPrm(2).Data;
n = robot.dof;

validateattributes(cfg.sample_time, {'double'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(cfg.q_limits, {'double'}, {'size', [n, 2], 'real', 'finite'});
assert(all(cfg.q_limits(:, 1) <= cfg.q_limits(:, 2)), ...
    'pose_ik_feedback_step:InvalidJointLimits', 'Lower limits must not exceed upper limits.');
validateattributes(cfg.damping, {'double'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(cfg.step_size, {'double'}, {'scalar', 'real', 'positive', '<=', 1});
validateattributes(cfg.tolerance, {'double'}, {'vector', 'numel', 2, 'real', 'finite', 'positive'});
if isfield(cfg, 'qd_max')
    validateattributes(cfg.qd_max, {'double'}, {'vector', 'real', 'nonnan', 'nonnegative', 'nonempty'});
    assert(isscalar(cfg.qd_max) || numel(cfg.qd_max) == n, ...
        'pose_ik_feedback_step:InvalidSpeedLimits', 'Specify one speed limit or one per joint.');
end

block.NumInputPorts = 2;
block.NumOutputPorts = 2;
block.InputPort(1).Dimensions = [4 4];
block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = 'Real';
block.InputPort(1).DirectFeedthrough = true;
block.InputPort(2).Dimensions = n;
block.InputPort(2).DatatypeID = 0;
block.InputPort(2).Complexity = 'Real';
block.InputPort(2).DirectFeedthrough = true;
block.OutputPort(1).Dimensions = n;
block.OutputPort(1).DatatypeID = 0;
block.OutputPort(1).Complexity = 'Real';
block.OutputPort(2).Dimensions = 1;
block.OutputPort(2).DatatypeID = 8;  % Boolean
block.OutputPort(2).Complexity = 'Real';
block.SampleTimes = [cfg.sample_time, 0];

block.RegBlockMethod('Outputs', @compute_outputs);
end

function compute_outputs(block)
robot = block.DialogPrm(1).Data;
cfg = block.DialogPrm(2).Data;
T_d = block.InputPort(1).Data;
q_measured = block.InputPort(2).Data(:);
assert(all(isfinite(q_measured)), 'pose_ik_feedback_step:InvalidFeedback', ...
    'Measured joint positions must be finite.');
[q_command, converged] = position_ik_step(robot, T_d, q_measured, cfg);
block.OutputPort(1).Data = q_command;
block.OutputPort(2).Data = converged;
end

function [q_command, converged] = position_ik_step(robot, T_d, q_measured, cfg)
q_command = q_measured;
converged = false;
R_d = T_d(1:3, 1:3);
if ~all(isfinite(T_d(:))) || norm(T_d(4, :) - [0 0 0 1], inf) > 1e-8 || ...
        norm(R_d.' * R_d - eye(3), 'fro') > 1e-6 || det(R_d) <= 0
    return;
end

% Lesson 6: use the SE(3) logarithm and Jacobian in the same body frame.
[J_b, T] = jacobian_matrix(robot, q_measured);
e_b = logT(tform_inv(T) * T_d).';
if ~all(isfinite(e_b))
    return;
end
converged = norm(e_b(1:3)) < cfg.tolerance(1) && ...
    norm(e_b(4:6)) < cfg.tolerance(2);
if converged
    return;
end
regularization = cfg.damping^2 * eye(6);
delta_q = J_b.' * ((J_b * J_b.' + regularization) \ e_b);
if all(isfinite(delta_q))
    q_next = q_measured + cfg.step_size * delta_q;
    q_command = clip(q_next, cfg.q_limits(:, 1), cfg.q_limits(:, 2));
    if isfield(cfg, 'qd_max')
        qd_max = cfg.qd_max(:);
        max_step = cfg.sample_time * qd_max;
        % Limit recovery steps even when feedback is outside the position limits.
        q_command = clip(q_command, q_measured - max_step, q_measured + max_step);
    end
end
end
