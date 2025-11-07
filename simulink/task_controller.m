function task_controller(block)
setup(block);
end

function setup(block)
% SETUP 函数定义了模块的基本特性。
block.NumDialogPrms = 3; % robot
% block.DialogPrmsTunable = {'Tunable', 'Tunable'};

%% ① 注册输入/输出端口的数量和属性
block.NumInputPorts  = 3; % 输入端口数量
block.NumOutputPorts = 1; % 输出端口数量

% 配置输入端口属性
% q
block.InputPort(1).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(1).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(1).DirectFeedthrough = true;   % 直接馈通 (重要！)

% Td
block.InputPort(2).Dimensions        = [4, 4];      % 端口维度 (1 为标量)
block.InputPort(2).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(2).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(2).DirectFeedthrough = true;   % 直接馈通 (重要！)

% Vd
block.InputPort(3).Dimensions        = 6;      % 端口维度 (1 为标量)
block.InputPort(3).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(3).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(3).DirectFeedthrough = true;   % 直接馈通 (重要！)

% 配置输出端口属性 qd_ref
block.OutputPort(1).DatatypeID       = 0;
block.OutputPort(1).Complexity       = 'Real';

%% ② 设置模块的采样时间
block.SampleTimes = [-1, 0];
% [-1, 0] 表示继承采样时间
% [0, 0]  表示连续采样时间
% [T, offset] 表示离散采样时间，周期为 T，偏移为 offset

block.NumContStates = 6;
block.SimStateCompliance = 'DefaultSimState';

%% ④ 注册模块将使用的回调方法
% 这是 Level-2 S-Function 的核心
% block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup); % 初始化 DWork 向量
block.RegBlockMethod('InitializeConditions', @InitConditions);  % 初始化条件
% block.RegBlockMethod('Start', @Start);                         % 模块启动时调用
block.RegBlockMethod('Outputs', @Outputs);                     % 计算输出
% block.RegBlock('Update', @Update);                             % 更新离散状态
block.RegBlockMethod('Derivatives', @Derivatives);             % 计算导数 (连续状态)
% block.RegBlockMethod('Terminate', @Terminate);                 % 仿真结束时调用
block.RegBlockMethod('SetInputPortDimensions',  @SetInputPortDims);
% block.RegBlockMethod('SetOutputPortDimensions', @SetOutputPortDims);

end

function SetInputPortDims(block, idx, di)
    block.InputPort(idx).Dimensions = di;
    % 假设输出与输入维度一致
    if idx == 1
        block.OutputPort(1).Dimensions = di;
    end
end

% function SetOutputPortDims(block, idx, di)
%     block.OutputPort(idx).Dimensions = di;
% end

% function DoPostPropSetup(block)
%    block.OutputPort(1).Dimensions = block.InputPort(1).Dimensions;
% end

function InitConditions(block)
% 初始化状态
block.ContStates.Data = zeros(6,1);
end

function Outputs(block)
robot = block.DialogPrm(1).Data;
Kp = block.DialogPrm(2).Data(:);
Ki = block.DialogPrm(3).Data(:);
q = block.InputPort(1).Data(:);
Td = block.InputPort(2).Data;
T = forward_kin_general(robot, q);
R = T(1:3, 1:3);
p = T(1:3, 4);
Rd = Td(1:3, 1:3);
pd = Td(1:3, 4);
xe = [logR(R'*Rd)'; R' * (pd - p)];
sum_xe = block.ContStates.Data;
Vd = blkdiag(R', R') *  block.InputPort(3).Data(:);
Jb = jacobian_matrix(robot, q);
V_ref = Vd + Kp .* xe + Ki .* sum_xe;
out = lsqminnorm(Jb,V_ref);
block.OutputPort(1).Data = out(:);
end

%%
function Derivatives(block)
robot = block.DialogPrm(1).Data;
q = block.InputPort(1).Data(:);
Td = block.InputPort(2).Data;
T = forward_kin_general(robot, q);
R = T(1:3, 1:3);
p = T(1:3, 4);
Rd = Td(1:3, 1:3);
pd = Td(1:3, 4);
block.Derivatives.Data = [logR(R'*Rd)'; R' * (pd - p)];
end