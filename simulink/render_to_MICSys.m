function render_to_MICSys(block)
% S-Function 用于通过 UDP 发送关节数据到 MICSys 系统
setup(block);
end

function setup(block)
    % 输入输出端口配置
    block.NumInputPorts  = 1;
    block.NumOutputPorts = 0;
    
    % 设置输入端口属性
    block.SetPreCompInpPortInfoToDynamic;
    block.InputPort(1).Dimensions        = -1;  % 继承维度
    block.InputPort(1).DirectFeedthrough = true;
    block.InputPort(1).DatatypeID        = 0;   % double
    block.InputPort(1).Complexity        = 'Real';
    
    % 采样时间配置 [-1 0] 表示继承采样时间
    block.SampleTimes = [-1 0];
    
    % 仿真状态兼容性
    block.SimStateCompliance = 'DefaultSimState';
    port = udpport("datagram");
    % 注册回调方法
    % block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
    % block.RegBlockMethod('InitializeConditions', @InitializeConditions);
    block.RegBlockMethod('Outputs', @(block) Output(block, port));
    block.RegBlockMethod('Terminate', @(block) Terminate(block, port));
    % block.RegBlockMethod('Start', @Start);
end



function Output(block, port)
    if block.IsMajorTimeStep
    jt =  block.InputPort(1).Data;
    set_joints(port, jt);
    end
end

function Terminate(block, port)
    clear port; 
end
