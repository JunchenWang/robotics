function simulate_DO_control

port = udpport("byte");
robot = convert_robot_tree2(importrobot('urdf\iiwa7\iiwa7.urdf'));
robot2 = robot;
robot2.mass = 1.2*robot.mass;% error
n = robot.dof;

Kx = 20 * eye(6);
Bx = 20 * eye(6);
Bn = diag(ones(1,n) * 5);
Kn = diag(ones(1,n) * 5);
tspan = [0, 10];
MassMatrix = @(t, y) [eye(n), zeros(n, 2 * n); zeros(n), mass_matrix(robot, y(1:n)), zeros(n); zeros(n, 2 * n), eye(n)];
opts = odeset('Mass',MassMatrix,'OutputFcn',@(t, y, flag) odeplot(t, y, flag, port, robot));

y0 = zeros(3*n,1);
y0(1:n) = [0 75 0 -94 0 -81 0] / 180 * pi;
ptp(port, y0(1:n)');
kesai = cal_kuka_kesai(y0);
Ts = forward_kin_general(robot, y0);
Te = Ts;
Te(3,4) = Te(3,4) + 0.5;

freq = 500;
N = tspan(2) * freq + 1;
[~,~,~,~,pp] = trapveltraj([0, 1], N, 'EndTime', tspan(2));
pp = pp{1};
p = @(t) desired_task_pos(t, Ts, Te, pp, fnder(pp, 1), fnder(pp, 2));

Y = 10 * eye(n);
controller = @(t, y) DO_controller(robot2, p, Kx, Bx, Bn, Kn, kesai, Y, t, y);
control_target = @(t, y) manipulator_dynamics_observer(robot, controller, @(t, y) Wrench(t, y, robot), t, y);
[t,y] = ode15s(control_target,tspan,y0,opts);
cnt = length(t);
torque = zeros(cnt, n);
pos_error = zeros(cnt, 1);
rot_error = zeros(cnt, 1); 
for i = 1 : cnt
    Xd = p(t(i));
    X = forward_kin_general(robot, y(i, 1:n)) ;
    pos_error(i) = norm(Xd(1:3,4) - X(1:3,4));
    rot_error(i) = norm(logR(X(1:3,1:3)' * Xd(1:3,1:3)));
    torque(i,:) = controller(t(i), y(i,:)');
end
disp(pos_error(end));
figure;
plot(t, y(:,2 * n + 1:end),'-', 'LineWidth', 2); % disturbance
xlabel("$t$/s", 'interpreter','latex');
ylabel('$\tau_d$/Nm', 'interpreter','latex');
% yticks([0,.2, .4, .6, .8, 1.0, 1.2, 1.4]);
xticks([0,1,2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
lg = legend('关节1','关节2','关节3','关节4','关节5','关节6','关节7','Orientation','horizontal');
fontsize(lg,18,'points')
set(gcf,'Position',[100 100 1200 800]);
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

function [tao, td] = DO_controller(robot, desired_pos, Kx, Bx, Bn, Kn, kesai, Y, t, y)
% Xd is desired motion in task space
n = robot.dof;
q = y(1:n);% + 1e-3 * rand(1, 7); % noise
qd = y(n + 1: 2 * n);% + 1e-3 * rand(1, 7); %noise
[M, C, G, Jb, dJb, dM, dX, X] = m_c_g_matrix(robot,q,qd);
R = X(1:3,1:3);
p = X(1:3,4);
Vb = Jb * qd;
wb = Vb(1:3);
v = R * Vb(4:6);
Jh = [eye(3), zeros(3); zeros(3), R];
dJh = [zeros(3), zeros(3); zeros(3), dX(1:3,1:3)];
dJb = dJh * Jb + Jh * dJb;
Jb = Jh * Jb;
td = y(2 * n + 1 : end);

[Xd, vel, acc] = desired_pos(t);
Rd = Xd(1:3, 1:3);
pd = Xd(1:3, 4);
wd = R' * vel(1:3);
vd = vel(4:6);
alphad = R' * acc(1:3);
ad = acc(4:6);

q0 = inverse_kin_kuka_robot_kesai_near(robot, Xd, kesai, q)';
if isempty(q0)
    error('no inverse');
end
Vd = [wd;vd];
dVd = [alphad - cross(wb, wd) ;ad];
xe = [logR(R'*Rd)'; pd - p];
dxe = Vd - [wb;v];

Z = null_z(Jb);
dZ = derivative_null_z(Jb, dJb);


% ax1 = dVd + Kx * xe + Bx * dxe;
% ax1 = dVd + A_x_inv(Jb, M) * ((Mu_x(Jb, M, dJb, C) + Bx) * dxe + Kx * xe); % PD+
s = dxe + Kx * xe;
ax1 = dVd + Kx * dxe + A_x_inv(Jb, M) * ((Mu_x(Jb, M, dJb, C) + Bx) * s );
a1 = pinv_J_x(Jb, M, ax1 - dJb * qd);

qe = q0 - q;
qed = -qd;
a2 = null_proj(Jb, M, M \ (Bn * qed + Kn * qe));
% ax2 = A_v(Z, M) \ ((Mu_v(Z, M, dZ, dM, C) + Bn(1)) * (-pinv_Z(Z, M) * qd) + Z' * Kn(1) * qe);
% a2 = Z * (ax2 - derivative_pinv_Z(Z, M, dZ, dM) * qd);

% tao = M * (a1 + a2) + C * qd + G;
% td = -Y * pinv_J_x(Jb, M, s);
P = Y * qd;
tao_d = td + P;
% tao_d = tao_d - M * Z * pinv_Z(Z, M) * tao_d;
tao_d = Jb' * ((Jb * (M \ Jb'))  \ (Jb * (M \ tao_d))); % minus the component projected into the null space !
% disp(tao_d);
tao = M * (a1 + a2) + C * qd + G - tao_d;
td = Y * (M \ (C * qd + G - tao - P - td));
end


function yd = manipulator_dynamics_observer(robot, controller, fext, t, y)
% fext is a force function applied to the robot
% fext(t,y)(:,end) is applied to TCP and others are applied to link frames
n = robot.dof;
yd = zeros(3*n,1);
q = y(1:n);
qd = y(n + 1 : 2 * n);
hqqd = gravity_velocity_torque(robot, q, qd);
yd(1:n) = qd;
ext_torque = get_ext_torque(robot, q, fext(t, y));
[tau, td] = controller(t, y);
yd(n + 1 : 2 * n) = tau - hqqd + ext_torque;
yd(2 * n + 1 : end) = td;
end


function F = Wrench(t, y, robot)
F = zeros(6, robot.dof);
if t > 1
    F(:,4) = [0, 0, 0, 0, 0, 10]';
    F(:,7) = [0, 0, 0, 0, 10, 0]';
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
