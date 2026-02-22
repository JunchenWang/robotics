function general_integrator(block)
setup(block);
end

function setup(block)
%% 一个参数，表示初始条件
block.NumDialogPrms = 1;
%% 注册输入/输出端口的数量和属性
block.NumInputPorts  = 1; % 输入端口数量
block.NumOutputPorts = 1; % 输出端口数量
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
% 配置输入端口属性
block.InputPort(1).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.InputPort(1).Complexity        = 'Real'; % 数值类型 (实数/复数)
block.InputPort(1).DirectFeedthrough = false;   % 直接馈通 
% 配置输出端口属性 
block.OutputPort(1).DatatypeID        = 0;      % 数据类型 (0 表示 double)
block.OutputPort(1).Complexity        = 'Real'; % 数值类型 (实数/复数)
%% 设置模块的采样时间
block.SampleTimes = [-1, 0]; % [-1, 0] 表示继承采样时间
block.NumContStates = 1; % tell model it has states
%% 注册模块将使用的回调方法
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup); % 初始化 DWork 向量
block.RegBlockMethod('InitializeConditions', @InitConditions);  % 初始化条件
block.RegBlockMethod('Outputs', @Outputs);                     % 计算输出
block.RegBlockMethod('Derivatives', @Derivatives);             % 计算导数 (连续状态)
block.RegBlockMethod('SetInputPortDimensions',  @SetInputPortDims);
block.RegBlockMethod('SetOutputPortDimensions', @SetOutputPortDims);
end

function SetInputPortDims(block, idx, di)
    block.InputPort(idx).Dimensions = di;
    block.OutputPort(idx).Dimensions = di;
end

function SetOutputPortDims(block, idx, di)
    block.InputPort(idx).Dimensions = di;
    block.OutputPort(idx).Dimensions = di;
end

function DoPostPropSetup(block)
   block.OutputPort(1).Dimensions = block.InputPort(1).Dimensions;
   in_shape = block.InputPort(1).Dimensions;
   block.NumContStates = prod(in_shape);
   block.NumDworks = 1;
   block.Dwork(1).Name            = 'shape';
   block.Dwork(1).DatatypeID      = 0;      % double
   block.Dwork(1).Complexity      = 'Real';
   block.Dwork(1).UsedAsDiscState = false;
   block.Dwork(1).Dimensions = numel(in_shape);
end

function InitConditions(block)
x0 = block.DialogPrm(1).Data(:);
in_shape = block.InputPort(1).Dimensions;
n = prod(in_shape);
if isscalar(x0)
    x0 = repmat(x0, n, 1);
end
block.ContStates.Data = x0;
block.Dwork(1).Data = in_shape;
end

function Outputs(block)
shape = block.Dwork(1).Data;
if isscalar(shape)
    block.OutputPort(1).Data = block.ContStates.Data;
else
    block.OutputPort(1).Data = reshape(block.ContStates.Data,  shape(:)');
end
end

function Derivatives(block)
u = block.InputPort(1).Data(:);
block.Derivatives.Data = u;
end