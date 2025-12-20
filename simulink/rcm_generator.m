function rcm_generator(block)
setup(block);
end

function setup(block)
% SETUP 函数定义了模块的基本特性。
block.NumDialogPrms = 5; % robot, Ts, p1, p2, lambda0
% block.DialogPrmsTunable = {'Tunable', 'Tunable'};
robot = block.DialogPrm(1).Data;
n = robot.dof;
%% ① 注册输入/输出端口的数量和属性
block.NumInputPorts  = 2; % 输入端口数量
block.NumOutputPorts = 3; % 输出端口数量

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% q
block.InputPort(1).Dimensions       = n;
block.InputPort(1).DatatypeID       = 0;
block.InputPort(1).Complexity       = 'Real';
block.InputPort(1).DirectFeedthrough = true;   % 直接馈通 (重要！)
% qd
block.InputPort(2).Dimensions       = n;
block.InputPort(2).DatatypeID       = 0;
block.InputPort(2).Complexity       = 'Real';
block.InputPort(2).DirectFeedthrough = true;   % 直接馈通 (重要！)

% 配置输出端口属性 xe
block.OutputPort(1).Dimensions = 6;
block.OutputPort(1).DatatypeID       = 0;
block.OutputPort(1).Complexity       = 'Real';
% 配置输出端口属性 xd_d
block.OutputPort(2).Dimensions = 6;
block.OutputPort(2).DatatypeID       = 0;
block.OutputPort(2).Complexity       = 'Real';
% 配置输出端口属性 J
block.OutputPort(3).Dimensions = [6, n];
block.OutputPort(3).DatatypeID       = 0;
block.OutputPort(3).Complexity       = 'Real';
%% ② 设置模块的采样时间
block.SampleTimes = [-1, 0];
% [-1, 0] 表示继承采样时间
% [0, 0]  表示连续采样时间
% [T, offset] 表示离散采样时间，周期为 T，偏移为 offset

% block.NumContStates = 6;
% block.SimStateCompliance = 'DefaultSimState';

%% ④ 注册模块将使用的回调方法
% 这是 Level-2 S-Function 的核心
% block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup); % 初始化 DWork 向量
% block.RegBlockMethod('InitializeConditions', @InitConditions);  % 初始化条件
% block.RegBlockMethod('Start', @Start);                         % 模块启动时调用
block.RegBlockMethod('Outputs', @Outputs);                     % 计算输出
% block.RegBlock('Update', @Update);                             % 更新离散状态
% block.RegBlockMethod('Derivatives', @Derivatives);             % 计算导数 (连续状态)
% block.RegBlockMethod('Terminate', @Terminate);                 % 仿真结束时调用
% block.RegBlockMethod('SetInputPortDimensions',  @SetInputPortDims);
% block.RegBlockMethod('SetOutputPortDimensions', @SetOutputPortDims);

end

% function SetInputPortDims(block, idx, di)
%     block.InputPort(idx).Dimensions = di;
%     % 假设输出与输入维度一致
%     if idx == 3
%         block.OutputPort(1).Dimensions = di(2);
%     end
% end

function Outputs(block)
%相对于初始位置 
robot = block.DialogPrm(1).Data;
t = block.CurrentTime;
Ts = block.DialogPrm(2).Data;
Rs = Ts(1:3,1:3);
ps = Ts(1:3,4);
p1 = block.DialogPrm(3).Data(:);
p2 = block.DialogPrm(4).Data(:);
lambda0 = block.DialogPrm(5).Data;
% xyz = [0;0.02*t; p2];
% xyz = [0.02*sin(2*t);0.02*t; p2];
xyz_d = ps + Rs * [0.05*cos(t); 0.05 * sin(t); p2];
xyzd_d = Rs * [-0.05*sin( t); 0.05 *cos(t); 0];
% xyz = [0.02*cos(t); 0.02 * sin(t); p2 + 0.005*t];
% xyzd = [-0.02*sin(t); 0.02 * cos(t); 0];
% xyzd = zeros(3,1);

P1 = ps + Rs * [0, 0, p1]';
P2 = ps + Rs * [0, 0, p2]';
prcm_d = P1 + (P2 - P1) * lambda0;% + [0,0,0.01]'; % change rcm


q = block.InputPort(1).Data(:);
qd = block.InputPort(2).Data(:);

[J, ~, x, err, xyz, J2] = rcm_jacobian(robot, q, qd, [0, 0, p1]', [0, 0, p2]', prcm_d);
% disp(err);
xe = [xyz_d - xyz; prcm_d - x];
block.OutputPort(1).Data = xe;
block.OutputPort(2).Data = [xyzd_d; 0;0;0];
block.OutputPort(3).Data = [J2; J];


end

