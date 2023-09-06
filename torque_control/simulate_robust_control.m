function simulate_robust_control

port = udpport("byte");
robot = convert_robot_tree2(importrobot('urdf\iiwa7\iiwa7.urdf'));
robot2 = robot;
robot2.mass = 1.2 * robot.mass;% error
param1 = get_dynamics(robot);
param2 = get_dynamics(robot2);
rho = 1.2 * norm(param1 - param2);
n = robot.dof;

tspan = [0, 10];
MassMatrix = @(t, y) [eye(n), zeros(n); 
                      zeros(n), mass_matrix(robot, y(1:n))];
opts = odeset('Mass',MassMatrix,'OutputFcn',@(t, y, flag) odeplot(t, y, flag, port));

y0 = zeros(2*n,1);
y0(1:n) = [0 75 0 -94 0 -81 0] / 180 * pi;
qs = y0(1:n);
ptp(port, qs');
kesai = cal_kuka_kesai(qs);
Ts = forward_kin_general(robot, qs);
Te = Ts;
Te(3,4) = Te(3,4) + 0.5;
qe = inverse_kin_kuka_robot_kesai_near(robot, Te, kesai, qs)';
freq = 500;
N = tspan(2) * freq + 1;
[~,~,~,~,pp] = trapveltraj([0, 1], N, 'EndTime', tspan(2));
pp = pp{1};
p = @(t) desired_joint_pos(t, qs, qs, pp, fnder(pp, 1), fnder(pp, 2));

% Y = 1 * eye(10 *n);
A = 20 * eye(n);
K = 20 * eye(n);
controller = @(t, y) robust_controller(robot2, p, rho, A, K, t, y);
control_target = @(t, y) manipulator_dynamics(robot, controller, @Wrench, t, y);
[t,y] = ode15s(control_target,tspan,y0,opts);
cnt = length(t);
torque = zeros(cnt, n);
error = zeros(cnt, n);
for i = 1 : cnt
    Xd = p(t(i))';
    X =  y(i, 1:n);
    error(i,:) = X - Xd;
    torque(i,:) = controller(t(i), y(i,:)');
end
disp(error(end,:));
% figure;
% plot(t, y(:,2 * n + 1 : end)); % disturbance
figure;
plot(t, torque);
figure;
plot(t, error);
figure;
plot(t, y(:,n + 1: 2 * n)); % speed


end

function ret = odeplot(t, y, flag, port)
if strcmp(flag, 'init') == 1
elseif isempty(flag)
    n = size(y,1) / 2;
    setJoints(port, y(1:n, end));
else
end
ret = 0;
end


function [q, qd, qdd] = desired_joint_pos(t, qs, qe, pp, ppd, ppdd)
% line with 
s = ppval(pp, t);
sd = ppval(ppd, t);
sdd = ppval(ppdd, t);
q = qs + (qe - qs)*s;
qd = (qe - qs)*sd;
qdd = (qe - qs)*sdd;
end

function [tao, dp] = robust_controller(robot, desired_pos, rho, A, K, t, y)
% Xd is desired motion in task space
n = robot.dof;
q = y(1:n);% + 1e-3 * rand(1, 7); % noise
qd = y(n + 1: 2 * n);% + 1e-3 * rand(1, 7); %noise
[dq, dqd, dqdd] = desired_pos(t);
qe = dq - q;
qed = dqd - qd;
r = qed + A * qe;
v = dqd + A * qe;
a = dqdd + A * qed;
[Y, YTr, M, C, G, Jb, dJb, dM, dX, X] = regressor_m_c_g_matrix(robot,q,qd,a,v,r);
eps = 0.1;
if norm(YTr) > eps
    dp = rho * YTr / norm(YTr);  
else
    dp = rho / eps * YTr;
end
param = get_dynamics(robot);
   
tao = Y* (param + dp) + K * r;
end


function yd = manipulator_dynamics(robot, controller, fext, t, y)
% fext is a force function applied to the robot
% fext(t,y)(:,end) is applied to TCP and others are applied to link frames
n = robot.dof;
yd = zeros(2*n,1);
q = y(1:n);
qd = y(n + 1 : 2 * n);
hqqd = gravity_velocity_torque(robot, q, qd);
yd(1:n) = qd;
ext_torque = get_ext_torque(robot, q, fext(t, y));
[tau, ~] = controller(t, y);
yd(n + 1 : 2 * n) = tau - hqqd + ext_torque;
end


function F = Wrench(t, y)
n = size(y,1) / 2;
F = zeros(6,n);
if t > 1
    % F(:,4) = [0, 0, 0, 0, 0, 10]';
end
end

function ptp(port, jts, vel)
if nargin < 3
    vel = .4;
end
start = queryJoints(port);
wayPoints = [start',jts'];
Freq = 200;
rate = rateControl(Freq);
T = max(abs(jts - start) / vel);
numSamples = round(T * Freq) + 1;
jt = trapveltraj(wayPoints,numSamples);
for i = 1 : numSamples
    setJoints(port, jt(:,i));
    waitfor(rate);
end
end

function setJoints(port, jt)
cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
    ,jt(6), jt(7));
writeline(port,cmd,"127.0.0.1",7755);
end

function joints = queryJoints(port)
% ";" 表示查询关节角
writeline(port,"robot;","127.0.0.1",7755);
s = readline(port);
joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
joints = mod(joints + pi, 2*pi) - pi;
end
