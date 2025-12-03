function robot_dynamics(block)
setup(block);
end

function setup(block)
% SETUP 函数定义了模块的基本特性。
block.NumDialogPrms = 3; % robot, init_pose, init_velocity
% block.DialogPrmsTunable = {'Tunable', 'Tunable'};
robot = block.DialogPrm(1).Data;
n = robot.dof;
%% ① 注册输入/输出端口的数量和属性
block.NumInputPorts  = 2; % 输入端口数量
block.NumOutputPorts = 2; % 输出端口数量

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% 配置输入端口属性
% tau
block.InputPort(1).Dimensions        = n;
block.InputPort(1).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(1).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(1).DirectFeedthrough = false;   % 直接馈通 (重要！)

% fext
block.InputPort(2).Dimensions        = [6, n];      % 端口维度 (1 为标量)
block.InputPort(2).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(2).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(2).DirectFeedthrough = false;   % 直接馈通 (重要！)


% 配置输出端口属性 q
block.OutputPort(1).Dimensions       = n;
block.OutputPort(1).DatatypeID       = 0;
block.OutputPort(1).Complexity       = 'Real';

% 配置输出端口属性 qd
block.OutputPort(2).Dimensions       = n;
block.OutputPort(2).DatatypeID       = 0;
block.OutputPort(2).Complexity       = 'Real';


%% ② 设置模块的采样时间
block.SampleTimes = [-1, 0];
% [-1, 0] 表示继承采样时间
% [0, 0]  表示连续采样时间
% [T, offset] 表示离散采样时间，周期为 T，偏移为 offset
block.NumContStates = 2 * n;
block.SimStateCompliance = 'DefaultSimState';

%% ④ 注册模块将使用的回调方法
% 这是 Level-2 S-Function 的核心
block.RegBlockMethod('CheckParameters', @CheckParam);
block.RegBlockMethod('InitializeConditions', @InitConditions);  % 初始化条件
block.RegBlockMethod('Outputs', @Outputs);                     % 计算输出
block.RegBlockMethod('Derivatives', @Derivatives);             % 计算导数 (连续状态)
% block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup); % 初始化 DWork 向量
% block.RegBlockMethod('Start', @Start);                         % 模块启动时调用
% block.RegBlock('Update', @Update);                             % 更新离散状态
% block.RegBlockMethod('Terminate', @Terminate);                 % 仿真结束时调用
% block.RegBlockMethod('SetInputPortDimensions',  @SetInputPortDims);
% block.RegBlockMethod('SetOutputPortDimensions', @SetOutputPortDims);

end


function CheckParam(block)
robot = block.DialogPrm(1).Data;
q = block.DialogPrm(2).Data;
qd = block.DialogPrm(3).Data;
if length(q) ~= robot.dof || length(qd) ~= robot.dof
    error('robot dof is not compatible with joint length');
end
end

function InitConditions(block)
% 初始化状态
init_q = block.DialogPrm(2).Data;
init_qd = block.DialogPrm(3).Data;
block.ContStates.Data = [init_q(:);init_qd(:)];
end

function Outputs(block)
robot = block.DialogPrm(1).Data;
n = robot.dof;
block.OutputPort(1).Data = block.ContStates.Data(1:n);
block.OutputPort(2).Data = block.ContStates.Data(n + 1: end);
end

function Derivatives(block)
robot = block.DialogPrm(1).Data;
n = robot.dof;
q = block.ContStates.Data(1:n);
qd = block.ContStates.Data(n + 1 : end);
tau = block.InputPort(1).Data(:);
fext = block.InputPort(2).Data;
% hqqd = gravity_velocity_torque(robot, q, qd);
% ext_torque = get_ext_torque(robot, q, fext);
M = mass_matrix(robot, q);
% tem = cqd + g - ext_torque
tem = inverse_dynamics_fext(robot, q, qd, zeros(n,1), fext);
block.Derivatives.Data = [qd; M \ (tau - tem)];
end

% function Start(block)
% robot = block.DialogPrm(1).Data;
% block.NumContStates = 2 * robot.dof;
% end

% function SetOutputPortDims(block, idx, di)
%     block.OutputPort(idx).Dimensions = di;
% end

% function DoPostPropSetup(block)
% block.NumDworks = 1;
% block.Dwork(1).Name = 'x0';
% block.Dwork(1).Dimensions      = 1;
% block.Dwork(1).DatatypeID      = 0;
% block.Dwork(1).Complexity      = 'Real';
% end


% function SetInputPortDims(block, idx, di)
% robot = block.DialogPrm(1).Data;
% if idx == 1 && prod(di) ~= robot.dof
%     error('tau input size is wrong');
% end
%
% if idx == 2 && (length(di) ~= 2 || di(2) ~= robot.dof)
%     error('external force input size is wrong');
% end
%
% block.InputPort(idx).Dimensions = di;
% if idx == 1
%     block.OutputPort(1).Dimensions = di;
%     block.OutputPort(2).Dimensions = di;
% end
% end