function pose_line_trajectory(block)
%POSE_LINE_TRAJECTORY Sampled straight-line reference with fixed orientation.
% Parameters: robot, world-frame displacement (m), duration (s), sample time (s).
% The first measured joint vector defines the trajectory's starting pose.

block.NumDialogPrms = 4;
block.DialogPrmsTunable = repmat({'Nontunable'}, 1, 4);
robot = block.DialogPrm(1).Data;
displacement = block.DialogPrm(2).Data;
duration = block.DialogPrm(3).Data;
sample_time = block.DialogPrm(4).Data;
validateattributes(displacement, {'double'}, {'vector', 'numel', 3, 'real', 'finite'});
validateattributes(duration, {'double'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(sample_time, {'double'}, {'scalar', 'real', 'finite', 'positive'});

block.NumInputPorts = 1;
block.NumOutputPorts = 1;
block.InputPort(1).Dimensions = robot.dof;
block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = 'Real';
block.InputPort(1).DirectFeedthrough = true;
block.OutputPort(1).Dimensions = [4 4];
block.OutputPort(1).DatatypeID = 0;
block.OutputPort(1).Complexity = 'Real';
block.SampleTimes = [sample_time, 0];

block.RegBlockMethod('PostPropagationSetup', @allocate_state);
block.RegBlockMethod('InitializeConditions', @initialize_state);
block.RegBlockMethod('Outputs', @compute_output);
block.RegBlockMethod('Update', @capture_start);
end

function allocate_state(block)
block.NumDworks = 2;
block.Dwork(1).Name = 'start_pose';
block.Dwork(1).Dimensions = 16;
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = 'Real';
block.Dwork(1).UsedAsDiscState = true;
block.Dwork(2).Name = 'start_time';
block.Dwork(2).Dimensions = 1;
block.Dwork(2).DatatypeID = 0;
block.Dwork(2).Complexity = 'Real';
block.Dwork(2).UsedAsDiscState = true;
end

function initialize_state(block)
block.Dwork(1).Data = zeros(16, 1);
block.Dwork(2).Data = NaN;
end

function compute_output(block)
start_time = block.Dwork(2).Data;
if isnan(start_time)
    block.OutputPort(1).Data = measured_pose(block);
    return;
end

displacement = block.DialogPrm(2).Data;
duration = block.DialogPrm(3).Data;
u = min(1, max(0, (block.CurrentTime - start_time) / duration));
s = u^3 * (10 + u * (-15 + 6*u));
T_d = reshape(block.Dwork(1).Data, 4, 4);
T_d(1:3, 4) = T_d(1:3, 4) + s * displacement(:);
block.OutputPort(1).Data = T_d;
end

function capture_start(block)
if isnan(block.Dwork(2).Data)
    T_start = measured_pose(block);
    block.Dwork(1).Data = T_start(:);
    block.Dwork(2).Data = block.CurrentTime;
end
end

function T = measured_pose(block)
robot = block.DialogPrm(1).Data;
q = block.InputPort(1).Data;
T = forward_kin_general(robot, q(:));
end
