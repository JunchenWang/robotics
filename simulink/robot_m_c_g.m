function robot_m_c_g(block)
setup(block);
end

function setup(block)
% SETUP 函数定义了模块的基本特性。
block.NumDialogPrms = 1; % robot, init_pose, init_velocity
% block.DialogPrmsTunable = {'Tunable', 'Tunable'};
robot = block.DialogPrm(1).Data;
n = robot.dof;
%% ① 注册输入/输出端口的数量和属性
block.NumInputPorts  = 2; % 输入端口数量
block.NumOutputPorts = 8; % 输出端口数量

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% 配置输入端口属性
% q
block.InputPort(1).Dimensions        = n;
block.InputPort(1).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(1).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(1).DirectFeedthrough = true;   % 直接馈通 (重要！)

% qd
block.InputPort(2).Dimensions        = n;      % 端口维度 (1 为标量)
block.InputPort(2).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(2).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(2).DirectFeedthrough = true;   % 直接馈通 (重要！)


% 配置输出端口属性 M
block.OutputPort(1).Dimensions       = [n, n];
block.OutputPort(1).DatatypeID       = 0;
block.OutputPort(1).Complexity       = 'Real';

% 配置输出端口属性 C
block.OutputPort(2).Dimensions       = [n, n];
block.OutputPort(2).DatatypeID       = 0;
block.OutputPort(2).Complexity       = 'Real';

% 配置输出端口属性 g
block.OutputPort(3).Dimensions       = n;
block.OutputPort(3).DatatypeID       = 0;
block.OutputPort(3).Complexity       = 'Real';

% 配置输出端口属性 Jb
block.OutputPort(4).Dimensions       = [6, n];
block.OutputPort(4).DatatypeID       = 0;
block.OutputPort(4).Complexity       = 'Real';

% 配置输出端口属性 dJb
block.OutputPort(5).Dimensions       = [6, n];
block.OutputPort(5).DatatypeID       = 0;
block.OutputPort(5).Complexity       = 'Real';

% 配置输出端口属性 dMq
block.OutputPort(6).Dimensions       = [n, n];
block.OutputPort(6).DatatypeID       = 0;
block.OutputPort(6).Complexity       = 'Real';

% 配置输出端口属性 dTcp
block.OutputPort(7).Dimensions       = [4, 4];
block.OutputPort(7).DatatypeID       = 0;
block.OutputPort(7).Complexity       = 'Real';

% 配置输出端口属性 Tcp
block.OutputPort(8).Dimensions       = [4, 4];
block.OutputPort(8).DatatypeID       = 0;
block.OutputPort(8).Complexity       = 'Real';
%% ② 设置模块的采样时间
block.SampleTimes = [-1, 0];
% [-1, 0] 表示继承采样时间
% [0, 0]  表示连续采样时间
% [T, offset] 表示离散采样时间，周期为 T，偏移为 offset
block.RegBlockMethod('Outputs', @Outputs); 
end


function Outputs(block)
robot = block.DialogPrm(1).Data;
q = block.InputPort(1).Data(:);
qd = block.InputPort(2).Data(:);
[Mq, C, g, Jb, dJb, dMq, dTcp, Tcp] = m_c_g_matrix(robot, q, qd);
block.OutputPort(1).Data = Mq;
block.OutputPort(2).Data = C;
block.OutputPort(3).Data = g;
block.OutputPort(4).Data = Jb;
block.OutputPort(5).Data = dJb;
block.OutputPort(6).Data = dMq;
block.OutputPort(7).Data = dTcp;
block.OutputPort(8).Data = Tcp;
end
