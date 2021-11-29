function kuka_control_simulation
u = udpport("byte");
robot = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
n = robot.dof;
if n ~= 7
    error('not kuka!');
end
pid = [2,0,4];
qd = @(t) [[pi/2 pi/3 pi/6 2*pi/3 -pi/2 -pi/3 0]'; zeros(n,1); zeros(n, 1)];
% qd = @(t) [[sin(t) 0 0 cos(t) 0 0 0]'; zeros(n,1); zeros(n, 1)];
tao = @(t, y) computed_torque_control(robot, qd, pid, t, y);
% tao2 = @(t, y) zeros(n,1);
MassMatrix = @(t, y) [eye(n), zeros(n); zeros(n), mass_matrix(robot, y(1:n))];
opts = odeset('Mass',MassMatrix,'OutputFcn',@odeplot);
y0 = zeros(2*n,1);
ptp(y0(1:n)');
clear computed_torque_control;
control_target = @(t, y) manipulator_dynamics(robot, tao, @ext_wrench,t, y); 
tspan = [0, 10];
[t,y] = ode45(control_target,tspan,y0,opts);
figure;
plot(t, y(:,1:n)');
    function Fext = ext_wrench(t, y)
        Fext = zeros(6,1);
    end
    function ret = odeplot(t, y, flag)
        if strcmp(flag, 'init') == 1
        elseif isempty(flag)
            setJoints(y(1:n));
        else
        end
        ret = 0;
    end
    function setJoints(jt)
        cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
            ,jt(6), jt(7));
        writeline(u,cmd,"192.168.3.34",7755);
    end
    function ptp(jts, vel)
        if nargin < 2
            vel = .4;
        end
        start = queryJoints;
        wayPoints = [start',jts'];
        Freq = 200;
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
        joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
        joints = mod(joints + pi, 2*pi) - pi;
    end
end
