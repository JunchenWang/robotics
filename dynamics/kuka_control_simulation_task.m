function kuka_control_simulation_task
u = udpport("byte");
robot = convert_robot_tree(importrobot('E:\data\URDF\iiwa7\iiwa7.urdf'));
robot2 = robot;
% robot2.mass(7) = 3.1;% error
n = robot.dof;
if n ~= 7
    error('not kuka!');
end
pid = [100,0,2, 5];


% refVel = zeros(6,1);
% refAcc = zeros(6,1);
% qd = @(t) [[sin(t) 0 0 cos(t) 0 0 0]'; zeros(n,1); zeros(n, 1)];
% tao = @(t, y) computed_torque_control_task(robot, refPose,refVel, refAcc, pid, t, y);
% tao2 = @(t, y) zeros(n,1);
MassMatrix = @(t, y) [eye(n), zeros(n); zeros(n), mass_matrix(robot, y(1:n))];
opts = odeset('Mass',MassMatrix,'OutputFcn',@odeplot);
y0 = zeros(2*n,1);
y0(1:7) = [0 75 0 -94 0 -81 0] / 180 * pi;
ptp(y0(1:n)');
T0 = forward_kin_kuka(y0);
T0(1:3,4) = T0(1:3,4) / 1000;
refPose = T0;
refPose(1,4) = refPose(1,4) + 0.15;
% refPose = trvec2tform([0.7 -.1 0.55]);
tInterval = [0, 10];
freq = 1000;
N = tInterval(2) * freq + 1;
tSamples = linspace(tInterval(1), tInterval(2), N);
[scale,sd,sdd] = trapveltraj([0, 1], N, 'EndTime', tInterval(2));
[tforms,v,a] = transformtraj(T0,refPose,tInterval,tSamples, 'TimeScaling', [scale;sd;sdd]);
fid = fopen('tform.txt', 'w');
fprintf(fid, '%.12f %.12f %.12f %.12f\n', tforms);
fclose(fid);
fid = fopen('v.txt', 'w');
fprintf(fid, '%f %f %f %f %f %f\n', v);
fclose(fid);
fid = fopen('a.txt', 'w');
fprintf(fid, '%f %f %f %f %f %f\n', a);
fclose(fid);
vel = zeros(6, N);
acc = zeros(6, N);
clear cartesian_velocity_estimator;

% fid = fopen('tform.txt', 'r');
% tforms = fscanf(fid, '%f %f %f %f', [4,inf]);
% fclose(fid);
% tforms = reshape(tforms, 4, 4,numel(tforms)/16);
for i = 1 : N
    [vel(:,i), acc(:,i)] = cartesian_velocity_estimator(tSamples(i), tforms(:,:,i));
end
tao = @(t, y) computed_torque_controller(robot2, tforms,vel, acc, pid, freq, t, y);
clear computed_torque_control_task;
control_target = @(t, y) manipulator_dynamics(robot, tao, @ext_wrench,t, y); 
torque = [];
[t,y] = ode23(control_target,tInterval,y0,opts);
figure;
plot(t, y(:,1:n)');
figure;
plot(torque');
    function Fext = ext_wrench(t, y)
        Fext = [0, 0, 0, 1, 0, 0]';
    end
    function ret = odeplot(t, y, flag)
        if strcmp(flag, 'init') == 1
        elseif isempty(flag)
            setJoints(y(1:n, end));
            torque = [torque, tao(t(end), y(:,end))];
        else
        end
        ret = 0;
    end
    function setJoints(jt)
        cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
            ,jt(6), jt(7));
        writeline(u,cmd,"127.0.0.1",7755);
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
        writeline(u,"robot;","127.0.0.1",7755);
        s = readline(u);
        joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
        joints = mod(joints + pi, 2*pi) - pi;
    end
end
