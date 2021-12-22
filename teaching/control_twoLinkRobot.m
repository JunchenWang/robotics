function control_twoLinkRobot
u = udpport("byte");
robot = twoLinkRobot;
%用于控制器的机器人参数，与理论值刻意造成偏差，因为实际动力学参数肯定是有误差的
robot2 = robot;
robot2.m = 2.35;
robot2.L = 0.55;
%机器人自由度
n = robot.dof;
%运动到初始位置
pos0 = [0, 0];
ptp(pos0);
% 机器人初始状态
y0 = [pos0, 0, 0]';
% 期望轨迹
q = @(t) [pi/2, pi/3]';
qd = @(t) [0;0];
qdd = @(t) [0;0];
% pid参数
kp = [800, 800];
kd = [130, 130];
ki = [0, 0];
%清除函数中的静态变量
clear torque_control;
%力矩控制器
tao = @(t, y) torque_control(robot2, q, qd, qdd, kp, ki, kd, t, y);
%机器人动力学模型
dynamic = @(t, y) twoLinkRobot_dynamics(t, y, robot, tao);
%仿真时间
tspan = [0, 10];
% ode函数的额外参数，每次解算都调用odeplot函数，传数据到仿真软件让机器人运动
opts = odeset('OutputFcn',@odeplot);
pret = 0;
%开始仿真
[t, y] = ode45(dynamic, tspan, y0, opts);

plot(t, y(:,1:2)');
    function ret = odeplot(t, y, flag)
        if strcmp(flag, 'init') == 1
        elseif isempty(flag)
            Freq = 1 / (t(end) - pret);
            pret = t(end);
            r = rateControl(Freq);
            waitfor(r);
            setJoints(y(1:n, end));
        else
        end
        ret = 0;
    end
    function setJoints(jt)
        cmd = sprintf('robot;%f;%f;', jt(1), jt(2));
        writeline(u,cmd,"192.168.3.34",7755);
    end

    function ptp(jts, vel)
        if nargin < 2
            vel = .4;
        end
        start = queryJoints;
        wayPoints = [start',jts'];
        Freq = 100;
        r = rateControl(Freq);
        T = max(abs(jts - start) / vel);
        numSamples = round(T * Freq) + 1;
        jt = trapveltraj(wayPoints,numSamples);
        for i = 1 : numSamples
            setJoints(jt(:,i));
            waitfor(r);
        end
    end
    function joints = queryJoints
        % ";" 表示查询关节角
        writeline(u,"robot;","192.168.3.34",7755);
        s = readline(u);
        joints = sscanf(s,'%f;%f;')';
        joints = mod(joints + pi, 2*pi) - pi;
    end

end
