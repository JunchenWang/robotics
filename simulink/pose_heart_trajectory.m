function pose_heart_trajectory(block)
%POSE_HEART_TRAJECTORY Sampled heart-shaped reference with fixed orientation.
% Parameters: robot, shape, duration (s), sample time (s).
% Shape can be a width (m) for an XZ heart starting at its lower tip, or a
% struct with width, plane ('xz' or 'yz'), and start_angle (rad).
% The first measured joint vector defines the trajectory's starting pose.

block.NumDialogPrms = 4;
block.DialogPrmsTunable = repmat({'Nontunable'}, 1, 4);
robot = block.DialogPrm(1).Data;
[width, ~, ~] = shape_parameters(block.DialogPrm(2).Data);
duration = block.DialogPrm(3).Data;
sample_time = block.DialogPrm(4).Data;
validateattributes(width, {'double'}, {'scalar', 'real', 'finite', 'positive'});
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

[width, axes, start_angle] = shape_parameters(block.DialogPrm(2).Data);
duration = block.DialogPrm(3).Data;
u = min(1, max(0, (block.CurrentTime - start_time) / duration));
s = u^3 * (10 + u * (-15 + 6*u));
T_d = reshape(block.Dwork(1).Data, 4, 4);
if u < 1
    offset = heart_coordinates(start_angle + 2*pi*s) - heart_coordinates(start_angle);
    T_d(axes, 4) = T_d(axes, 4) + (width/32)*offset;
end
block.OutputPort(1).Data = T_d;
end

function [width, axes, start_angle] = shape_parameters(shape)
width = shape;
axes = [1 3];
start_angle = pi;
if isstruct(shape)
    width = shape.width;
    plane = validatestring(shape.plane, {'xz', 'yz'});
    if strcmp(plane, 'yz')
        axes = [2 3];
    end
    start_angle = shape.start_angle;
    validateattributes(start_angle, {'double'}, {'scalar', 'real', 'finite'});
end
end

function point = heart_coordinates(theta)
point = [16*sin(theta)^3; ...
    13*cos(theta) - 5*cos(2*theta) - 2*cos(3*theta) - cos(4*theta)];
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
