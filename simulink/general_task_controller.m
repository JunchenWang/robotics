function general_task_controller(block)
setup(block);
end

function setup(block)
% SETUP 函数定义了模块的基本特性。
block.NumDialogPrms = 2; % Kp, Ki
block.DialogPrmsTunable = {'Tunable', 'Tunable'};

%% ① 注册输入/输出端口的数量和属性
block.NumInputPorts  = 3; % 输入端口数量
block.NumOutputPorts = 1; % 输出端口数量

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
% 配置输入端口属性
% xe
block.InputPort(1).Dimensions        = -1;
block.InputPort(1).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(1).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(1).DirectFeedthrough = true;   % 直接馈通 (重要！)

% xd_d
block.InputPort(2).Dimensions        = -1;
block.InputPort(2).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(2).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(2).DirectFeedthrough = true;   % 直接馈通 (重要！)

% J
block.InputPort(3).Dimensions        = [-1, -1];
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

% block.NumContStates = 6;
block.SimStateCompliance = 'DefaultSimState';

%% ④ 注册模块将使用的回调方法
% 这是 Level-2 S-Function 的核心
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup); % 初始化 DWork 向量
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
    if idx == 3
        block.OutputPort(1).Dimensions = di(2);
    end
end

% function SetOutputPortDims(block, idx, di)
%     block.OutputPort(idx).Dimensions = di;
% end

function DoPostPropSetup(block)
    block.NumContStates = block.InputPort(1).Dimensions(1);
   % block.OutputPort(1).Dimensions = block.InputPort(1).Dimensions;
end

function InitConditions(block)
% 初始化状态

block.ContStates.Data = zeros(block.NumContStates,1);
end

function Outputs(block)
Kp = block.DialogPrm(1).Data(:);
Ki = block.DialogPrm(2).Data(:);

xe = block.InputPort(1).Data(:);
xd_d = block.InputPort(2).Data(:);
J = block.InputPort(3).Data;
sum_xe = block.ContStates.Data;
xd = xd_d + Kp .* xe + Ki .* sum_xe;
qd = damping_least_square(J, xd, 0.05, 1e3);
block.OutputPort(1).Data = qd;
end

%%
function Derivatives(block)

xe = block.InputPort(1).Data(:);
block.Derivatives.Data = xe;
end