function simulate_dual_loop_admittance_control

port = udpport("byte");
robot = convert_robot_tree2(importrobot('urdf\iiwa7\iiwa7.urdf'));
robot2 = robot;
robot2.mass = robot.mass;% error
n = robot.dof;
Kx = zeros(6,6,3);
Bx = zeros(6,6,3);

task_choice = 2;% not change
null_choice = 3;
rep_choice = 2; % xe,dxe expression 1: couple 2: uncouple
Kx(:,:,1) = 10 * eye(6);% pd
Bx(:,:,1) = 10 * eye(6);

Kx(:,:,2) = 1000 * eye(6);% pd+, impedance
Bx(:,:,2) = 200 * eye(6);

Kx(:,:,3) = 5 * eye(6);% passivity
Bx(:,:,3) = 5 * eye(6);

Bn = 2 * eye(n);
Kn = 20 * eye(n);
tspan = [0, 10];
MassMatrix = @(t, y) [eye(n), zeros(n, 2 * n + 12); 
                        zeros(n), mass_matrix(robot, y(1:n)), zeros(n, n + 12); 
                        zeros(n, 2 * n), eye(n), zeros(n, 12);
                        zeros(12, 3 * n), eye(12)];
opts = odeset('Mass',MassMatrix,'OutputFcn',@(t, y, flag) odeplot(t, y, flag, port, robot));

y0 = zeros(3*n +12,1);
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
Y = 100 * eye(n);
M = [1,1,1,1,1,1]';
K = [1000,1000,1000,1000,1000,1000]';
B = 2*sqrt(K);
admittance = @(t, y, fext) addmittance_controller(robot2, motion_planner, M, B, K, t, y, fext);
controller = @(t, y, fext) DO_position_controller(robot2, admittance, rep_choice, task_choice, null_choice, Kx, Bx, Bn, Kn, kesai, Y, t, y, fext);
control_target = @(t, y) manipulator_dynamics_observer(robot, controller, @(t, y) Wrench(t, y, refZ, stiffness, robot), t, y);
[t,y] = ode15s(control_target,tspan,y0,opts);
cnt = length(t);
torque = zeros(cnt, n);
pos_error = zeros(cnt, 1);
rot_error = zeros(cnt, 1);
phi = zeros(cnt,1);
force = zeros(cnt, 6);
for i = 1 : cnt
    Xd = motion_planner.desired_pose(t(i),y(i,:)');
    X = forward_kin_general(robot, y(i, 1:n)) ;
    phi(i) = cal_kuka_kesai(y(i, 1:n));
    pos_error(i) = norm(Xd(1:3,4) - X(1:3,4));
    rot_error(i) = norm(logR(X(1:3,1:3)' * Xd(1:3,1:3)));
    
    W = Wrench(t(i), y(i,:)', refZ, stiffness, robot);
    torque(i,:) = controller(t(i), y(i,:)', W);
    X(1:3,4) = zeros(3,1);
    force(i,:) = -adjoint_T(tform_inv(X))'* W(:,end);% applied to env
end
disp(pos_error(end));

figure;
plot(t, y(:,2 * n + 1: 3 * n),'-', 'LineWidth', 2); % disturbance
xlabel("$t$/s", 'interpreter','latex');
ylabel('$\hat{\tau}_d$/Nm', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
xticks([0,1,2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
lg = legend('关节1','关节2','关节3','关节4','关节5','关节6','关节7','Orientation','horizontal');
fontsize(lg,18,'points')
set(gcf,'Position',[100 100 1200 800]);
savefig('干扰力矩.fig');
saveas(gcf, '干扰力矩.png');

figure;
plot(t, torque,'-', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('$\tau$/Nm', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
xticks([0,1,2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
lg = legend('关节1','关节2','关节3','关节4','关节5','关节6','关节7','Orientation','horizontal');
fontsize(lg,18,'points')
set(gcf,'Position',[100 100 1200 800]);
savefig('控制力矩.fig');
saveas(gcf, '控制力矩.png');

figure;
plot(t, rot_error,'-', 'LineWidth', 2);
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
plot(t, pos_error,'-', 'LineWidth', 2);
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
plot(t, phi,'-', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('臂角/rad', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
xticks([0,1,2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
set(gcf,'Position',[100 100 1200 800]);
savefig('臂角.fig');
saveas(gcf, '臂角.png');

figure;
plot(t, force(:,end),'-', 'LineWidth', 2);
xlabel("$t$/s", 'interpreter','latex');
ylabel('接触力/N', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
xticks([0,1,2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
set(gcf,'Position',[100 100 1200 800]);
savefig('接触力.fig');
saveas(gcf, '接触力.png');

end


function [Td, vel, acc] = desired_task_pos(t, Ts, Te, pp, ppd, ppdd)
% line with
s = ppval(pp, t);
sd = ppval(ppd, t);
sdd = ppval(ppdd, t);
ps = Ts(1:3,4);
pe = Te(1:3,4);
Rs = Ts(1:3,1:3);
Re = Te(1:3,1:3);
pd = ps + (pe - ps)*s;
xe = logR(Rs'*Re)';
Rd = Rs * exp_w(xe*s);
Td = [Rd, pd; 0 0 0 1];
vel = [Rs * xe; (pe - ps)] * sd;
acc = [Rs * xe; (pe - ps)] * sdd;
end

function [Td, vel, acc] = desired_task_pos_sin(t, Ts, Te, pp, ppd, ppdd)
% line with
s = ppval(pp, t);
sd = ppval(ppd, t);
sdd = ppval(ppdd, t);
ps = Ts(1:3,4);
pe = Te(1:3,4);
Rs = Ts(1:3,1:3);
Re = Te(1:3,1:3);
len = norm(pe - ps);%y方向距离
width = 0.05;%x方向振幅
y = len*s;
w = 50;%角速度
x = width * sin(w * y);
yd = len*sd;
ydd = len*sdd;
xd = width * w * cos(w * y) * yd;
xdd = -width * w^2 * sin(w * y) * yd^2 + width * w * cos(w * y) * ydd;
pd = ps + [x, y, 0]';
xe = logR(Rs'*Re)';
Rd = Rs * exp_w(xe*s);
Td = [Rd, pd; 0 0 0 1];
vel = [Rs * xe * sd; xd;yd;0];
acc = [Rs * xe * sdd; xdd;ydd;0];
end

function [Td, vel, acc, state] = addmittance_controller(robot, motion_planner, M, B, K, t, y, fext)
    n = robot.dof;
    Fex = -fext(:,end);
    state = zeros(12,1);

    [Td, vel, acc] = motion_planner.desired_pose(t, y);
    Rd = Td(1:3,1:3);
    pd = Td(1:3,4);
    wd = vel(1:3);
    vd = vel(4:6);
    alphad = acc(1:3);
    ad = acc(4:6);
    re = y(3 * n + 1 : 3 * n + 3);
    red = y(3 * n + 4 : 3 * n + 6);
    pe = y(3 * n + 7 : 3 * n + 9);
    ped = y(3 * n + 10 : 3 * n + 12);
    redd = (Fex(1:3) - K(1:3) .* re - B(1:3) .* red) ./ M(1:3);
    pedd = (Fex(4:6) - K(4:6) .* pe - B(4:6) .* ped) ./ M(4:6);
    state(1:3) = red;
    state(4:6) = redd;
    state(7:9) = ped;
    state(10:12) = pedd;
    
    R = Rd * exp_w(-re);
    p = pd - R * pe;
    A = w_dr_A(re);
    dA = derivative_Ar(re, red);
    w = wd - Rd * A * red;
    alpha = alphad - Rd * A * (redd + (A \ dA) * (A \ Rd') * (wd - w)) + cross(wd, w);
    v = vd - R * ped - cross(w, pd - p);
    a = ad - R * pedd - cross(alpha, pd - p) - 2 * cross(w, vd - v) + cross(w, cross(w, pd - p));
    vel = [w;v];
    acc = [alpha; a];

    % vel = zeros(6,1);
    % acc = zeros(6,1);
    Td = [R, p; 0 0 0 1];
end

function [tao, td] = DO_position_controller(robot, addmittance_controller, rep_choice, task_choice, null_choice, Kx, Bx, Bn, Kn, kesai, Y, t, y, fext)
% Xd is desired motion in task space
n = robot.dof;
q = y(1:n);% + 1e-3 * rand(1, 7); % noise
qd = y(n + 1: 2 * n);% + 1e-3 * rand(1, 7); %noise
td = y(2 * n + 1 : 3 * n);% tao_d disturbance
[M, C, G, Jb, dJb, dM, dX, X] = m_c_g_matrix(robot,q,qd);
[Xd, vel, acc, state] = addmittance_controller(t, y, fext);
R = X(1:3,1:3);
p = X(1:3,4);
Vb = Jb * qd;
wb = Vb(1:3);
v = R * Vb(4:6);
Rd = Xd(1:3, 1:3);
pd = Xd(1:3, 4);
%% 任务空间
if rep_choice == 1 % 耦合
    wd = Rd' * vel(1:3);
    vd = Rd' * vel(4:6);
    alphad = Rd' * acc(1:3);
    ad = Rd' * acc(4:6);
    Vd = [wd;vd];
    dVd = [alphad; ad - cross(wd, vd)];
    dXd = Xd * se_twist(Vd);
    [dXinv, Xinv] = derivative_tform_inv(X, dX);
    Xe = Xinv * Xd;
    dXe = dXinv * Xd + Xinv * dXd;
    [dAdXe, AdXe] = derivative_adjoint_T(Xe, dXe);
    xe = logT(Xe)';
    dxe = AdXe * Vd - Vb;
    dVd = dAdXe * Vd + AdXe * dVd;
else % 解耦
    Jh = [eye(3), zeros(3); zeros(3), R];
    dJh = [zeros(3), zeros(3); zeros(3), dX(1:3,1:3)];
    dJb = dJh * Jb + Jh * dJb;
    Jb = Jh * Jb;
    wd = R' * vel(1:3);
    vd = vel(4:6);
    alphad = R' * acc(1:3);
    ad = acc(4:6);
    Vd = [wd;vd];
    dVd = [alphad - cross(wb, wd) ;ad];
    xe = [logR(R'*Rd)'; pd - p];
    dxe = Vd - [wb;v];
end
%% 零空间
Z = null_z(Jb);
dZ = derivative_null_z(Jb, dJb);
if task_choice == 1 % pd
    ax1 = dVd + Kx(:,:,1) * xe + Bx(:,:,1) * dxe;
elseif task_choice == 2 % pd+
    ax1 = dVd + A_x_inv(Jb, M) * ((Mu_x(Jb, M, dJb, C) + Bx(:,:,2)) * dxe + Kx(:,:,2) * xe); % PD+
else % passivity
    s = dxe + Kx(:,:,3)  * xe;
    ax1 = dVd + Kx(:,:,3) * dxe + A_x_inv(Jb, M) * ((Mu_x(Jb, M, dJb, C) + Bx(:,:,3)) * s );
end
tao_x = M * pinv_J_x(Jb, M, ax1 - dJb * qd);

q0 = inverse_kin_kuka_robot_kesai_near(robot, Xd, kesai, q);
% if isempty(q0)
%     error('no inverse');
% end
qe = q0 - q;
qed = -qd;
if null_choice == 1
    a2 = null_proj(Jb, M, M \ (Bn * qed + Kn * qe));
    tao_n = M * a2;
elseif null_choice == 2
    a2 = null_proj(Jb, M, Bn * qed + Kn * qe);
    tao_n = M * a2;
elseif null_choice == 3
    ax2 = A_v(Z, M) \ ((Mu_v(Z, M, dZ, dM, C) + Bn(1)) * (-pinv_Z(Z, M) * qd) + Z' * Kn(1) * qe);
    a2 = Z * (ax2 - derivative_pinv_Z(Z, M, dZ, dM) * qd);
    tao_n = M * a2;
else
    tao0 = Bn * qed + Kn * qe;
    tao_n = tao0 - Jb' * ((Jb * (M \ Jb')) \ (Jb * (M \ tao0)));
end
%% 干扰力矩估计
P = Y * qd;
tao_d = td + P;
tao_d = Jb' * ((Jb * (M \ Jb'))  \ (Jb * (M \ tao_d))); % minus the component projected into the null space !
% disp(tao_d);
%% 求和
tao = tao_x + tao_n + C * qd + G - tao_d;
td = [Y * (M \ (C * qd + G - tao - P - td)); state];
end


function yd = manipulator_dynamics_observer(robot, controller, fext, t, y)
% fext is a force function applied to the robot
% fext(t,y)(:,end) is applied to TCP and others are applied to link frames
n = robot.dof;
yd = zeros(3*n + 12,1);
q = y(1:n);
qd = y(n + 1 : 2 * n);
hqqd = gravity_velocity_torque(robot, q, qd);
yd(1:n) = qd;
Fext = fext(t, y);
ext_torque = get_ext_torque(robot, q, Fext);
[tau, td] = controller(t, y, Fext);
yd(n + 1 : 2 * n) = tau - hqqd + ext_torque;
yd(2 * n + 1 : end) = td;
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
if t > 5
    % F(:,4) = [0,0,0,0,0,10]';
end
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
