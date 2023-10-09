function simulate_discrete_dual_loop_admittance_control_with_velocity

port = udpport("byte");
robot = convert_robot_tree2(importrobot('urdf\iiwa7\iiwa7.urdf'));
n = robot.dof;

Kp_s = [100,100,100,100,100,100,100]';
Ki_s = [50,50,50,50,100,100,100]';

Kp_p = 5 * [1, 1, 1, 1, 1, 1]';
Ki_p = 0 * [1,1,1,1,1,1]';
Kd_p = 0 * [1,1,1,1,1,1]';

d = @(t, y) [7;6;5;4;3;2;1]*10;
Jm = [8;5;4;3;2;2;2];
Bm = 2;
r = 200;

y0 = zeros(3 * n + 6, 1);

tspan = [0, 10];

opts = odeset('OutputFcn',@(t, y, flag) odeplot(t, y, flag, port, robot));


y0(1:n) = [-40 75 0 -94 0 -81 0] / 180 * pi;
Ts = forward_kin_general(robot, y0);
kesai = cal_kuka_kesai(y0);
Ts(1:3,1:3) = [-1, 0, 0; 0 1 0;0 0 -1];
y0(1:n) = inverse_kin_kuka_robot_kesai_near(robot, Ts, kesai, y0(1:n));
ptp(port, y0(1:n)');
Ts = forward_kin_general(robot, y0);
refZ = Ts(3,4);
Tds = cell(1,3);
Tds{1} = Ts;
Te = Ts;
Te(3,4) = Te(3,4) - 0.01;
Tds{2} = Te;
Te(2,4) = Te(2,4) + 0.8;
Tds{3} = Te;
tnodes = [0, 2, tspan(2)];
stiffness = 1e5; % N/m
motion_planner = LinearTrajectory(tnodes, Tds);


M = [1,1,1,1,1,1]';
K = [1000,1000,1000,1000,1000,1000]';
B = 2*sqrt(K);
dt = 0.002;
cnt = round(tspan(2) / dt) + 1;
re = zeros(3,1);
pe = zeros(3,1);
red = zeros(3,1);
ped = zeros(3,1);
y = y0;
pos_error = zeros(cnt, 1);
rot_error = zeros(cnt, 1);
phi = zeros(cnt,1);
force = zeros(cnt, 6);
tt = zeros(cnt, 1);
tic;
for i = 1 : cnt
    t = (i - 1) * dt;
    tt(i) = t;
    [Td, vel, acc] = motion_planner.desired_pose(t, y);
    %%前馈
    Rd = Td(1:3,1:3);
    pd = Td(1:3,4);
    wd = vel(1:3);
    vd = vel(4:6);
    alphad = acc(1:3);
    ad = acc(4:6);
    Vd = [Rd' * wd; Rd' * vd];
    Td_ = Td * exp_twist(Vd * dt);
    Rd_ = Td_(1:3,1:3);
    pd_ = Td_(1:3,4);
    wd_ = wd + alphad * dt;
    vd_ = vd + ad * dt;
    % X = forward_kin_general(robot2, y(1:n)) ;
    %% 反馈
    [J, X] = jacobian_matrix(robot, y(1:n));
    Vb = J * y(n + 1 : 2 * n);
    R = X(1:3,1:3);
    p = X(1:3,4);
    w = R * Vb(1:3);
    v = R * Vb(4:6);

    re = logR(R' * Td(1:3,1:3))';
    pe = R' * (Td(1:3,4) - p);
    red = w_dr_A(logR(R' * Td(1:3,1:3))') \ Td(1:3,1:3)' * (vel(1:3) - w);
    ped = R' * (vel(4:6) - v - cross(w, Td(1:3,4) - p));

    fext =  Wrench(t, y, refZ, stiffness, robot);
    Fex = -fext(:,end);
    redd = (Fex(1:3) - K(1:3).*re - B(1:3).*red) ./ M(1:3);
    re = re + red * dt + 0.5 * redd * dt^2;
    red = red + redd * dt;

    pedd = (Fex(4:6) - K(4:6).*pe - B(4:6).*ped) ./ M(4:6);
    pe = pe + ped * dt + 0.5 * pedd * dt^2;
    ped = ped + pedd * dt;


    R = Rd_ * exp_w(-re);
    p = pd_ - R * pe;
    A = w_dr_A(re);
    dA = derivative_Ar(re, red);
    w = wd_ - Rd_ * A * red;
    alpha = alphad - Rd_ * A * (redd + (A \ dA) * (A \ Rd_') * (wd_ - w)) + cross(wd_, w);
    v = vd_ - R * ped - cross(w, pd_ - p);
    a = ad - R * pedd - cross(alpha, pd_ - p) - 2 * cross(w, vd_ - v) + cross(w, cross(w, pd_ - p));
    vel = [w;v];
    % vel = zeros(6,1);
    acc = [alpha; a];
    Tr = [R, p; 0 0 0 1];




    phi(i) = cal_kuka_kesai(y(1:n));
    pos_error(i) = norm(Td(1:3,4) - X(1:3,4));
    rot_error(i) = norm(logR(X(1:3,1:3)' * Td(1:3,1:3)));
    X(1:3,4) = zeros(3,1);
    force(i,:) = -adjoint_T(tform_inv(X))'* fext(:,end);% applied to env
    v = @(t, y) task_pos_controller(t, y, Tr, vel, robot, Kp_p, Ki_p, Kd_p);
    u = @(t, y) speed_controller(t, y, v, Kp_s, Ki_s);
    dynamic = @(t, y) joint_motor_dynamic(t, y, u, d, Jm, Bm, r);

    [t, y] = ode45(dynamic, [t, t + dt], y, opts);
    y = y(end,:)';
end
toc;
disp(pos_error(end));



figure;
plot(tt, rot_error,'-', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('$||r_e||$/rad', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
xticks([0,1,2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
% lg = legend('关节1','关节2','关节3','关节4','关节5','关节6','关节7','Orientation','horizontal');
% fontsize(lg,18,'points')
set(gcf,'Position',[100 100 1200 800]);
savefig('姿态误差.fig');
saveas(gcf, '姿态误差.png');

figure;
plot(tt, pos_error,'-', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('$||p_e||$/m', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
xticks([0,1,2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
% lg = legend('关节1','关节2','关节3','关节4','关节5','关节6','关节7','Orientation','horizontal');
% fontsize(lg,18,'points')
set(gcf,'Position',[100 100 1200 800]);
savefig('位置误差.fig');
saveas(gcf, '位置误差.png');

figure;
plot(tt, phi,'-', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('臂角/rad', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
xticks([0,1,2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
set(gcf,'Position',[100 100 1200 800]);
savefig('臂角.fig');
saveas(gcf, '臂角.png');

figure;
plot(tt, force(:,end),'-', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('接触力/N', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
xticks([0,1,2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
set(gcf,'Position',[100 100 1200 800]);
savefig('接触力.fig');
saveas(gcf, '接触力.png');

end

function F = Wrench(t, y, zref, stiffness, robot)
F = zeros(6, robot.dof);
% F(:,4) = [0, 0, 0, 0, 0, 10]';
X = forward_kin_general(robot, y(1:robot.dof));
R = X(1:3,1:3);
z = R(:,3);
real_f = (X(3,4) - zref) * stiffness;
if real_f < 0
    X(1:3,4) = zeros(3,1); %增加约束力矩
    F(:,end) = adjoint_T(X)' * [100 * cross(z, [0,0,-1]');0;0;-real_f]; % from s to b
end
end


function [u, e] = speed_controller(t, y, desired_speed, Kp, Ki)
% y(1): theta, y(2) : w, y(3) sum e
n = (length(y) + 1) / 4;
[wd, e2] = desired_speed(t, y);
e1 = wd - y(n + 1: 2 * n);
u = Kp .* e1 + Ki .* y(2 * n + 1 : 3 * n);
e = [e1;e2];
end

function [desired_speed, e] = task_pos_controller(t, y, Td, vel, robot, Kp, Ki, Kd)
n = (length(y) + 1) / 4;
q = y(1:n);
qd = y(n + 1 : 2*n);
Jb = jacobian_matrix(robot, q);
V = Jb * qd;
T = forward_kin_general(robot, q);
R = T(1:3,1:3);
p = T(1:3,4);
wd = vel(1:3);
vd = vel(4:6);
Rd = Td(1:3,1:3);
pd = Td(1:3,4);
xe = logR(R'*Rd)';
pe = pd - p;
wb = R' * wd + Kp(1:3) .* xe + Ki(1:3) .* y(3*n + 1 : 3*n + 3) + Kd(1:3) .* (R' * wd - V(1:3));
vb = R'*(vd + Kp(4:6) .* pe + Ki(4:6) .* y(3*n + 4 : 3*n + 6) + Kd(4:6) .* (vd - R * V(4:6)));
e = [xe; pe];
desired_speed = lsqminnorm(Jb,[wb;vb]);
% disp(Jb * pinv(Jb));
end


function yd = joint_motor_dynamic(t, y, u, d, J, B, r)
% y(1): theta, y(2): w, y(3): sum e
n = (length(y) + 1) / 4;
yd = zeros(3 * n + 6, 1);
[U, e] = u(t,y);
yd(1:n) = y(n +1 : 2 * n);
yd(n +1 : 2*n) = (U - d(t, y) ./ r - B .* y(n +1: 2 *n)) ./ J;
yd(2 * n + 1 : end) = e;
% disp(u(t,y));
end


function ret = odeplot(t, y, flag, port, robot)
if strcmp(flag, 'init') == 1
elseif isempty(flag)
    setJoints(port, y(1:robot.dof, end));
else
end
ret = 0;
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
cmd = 'robot;' + join(string(jt),';') + ';';
writeline(port,cmd,"127.0.0.1",7755);
end

function joints = queryJoints(port)
writeline(port,"robot;","127.0.0.1",7755);
s = split(readline(port), ';');
joints = double(s(1:end-1))';
end
