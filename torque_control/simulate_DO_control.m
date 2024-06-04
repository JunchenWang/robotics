function simulate_DO_control

port = udpport("byte");
robot = convert_robot_tree2(importrobot('urdf\iiwa7\iiwa7.urdf'));
nominal_robot = robot;
nominal_robot.mass = 1.2*robot.mass;% error
n = robot.dof;
Kx = 3000 * eye(6);
Bx = 300 * eye(6);
Bn = 1 * eye(n);
Kn = 1 * eye(n);
tspan = [0, 10];
MassMatrix = @(t, y) [eye(n), zeros(n, 2 * n); zeros(n), mass_matrix(robot, y(1:n)), zeros(n); zeros(n, 2 * n), eye(n)];
opts = odeset('Mass',MassMatrix,'OutputFcn',@(t, y, flag) odeplot_micsys(t, y, flag, port, robot));

y0 = zeros(3*n,1);
y0(1:n) = [0 75 0 -94 0 -81 0] / 180 * pi;
ptp_move(port, y0(1:n)');

kesai = cal_kuka_kesai(y0);
Ts = forward_kin_general(robot, y0);
Te = Ts;
Te(3,4) = Te(3,4) + 0.5;
Td = {Ts, Te};
traj = LinearTrajectory(tspan, Td);
p = @(t) traj.desired_pose(t);

Y = 10 * eye(n);
DOB = @(robot, tau, M, C, G, J, y) manipulator_DOB(robot, Y, tau, M, C, G, J, y);
controller = @(t, y, Fext) PD_controller(nominal_robot, p, Kx, Bx, Kn, Bn, kesai, DOB, t, y);
target_sysm = @(t, y) manipulator_dynamics_general(robot, controller, @(t, y) Wrench(t, y, robot), t, y);
[t,y] = ode15s(target_sysm,tspan,y0,opts);
cnt = length(t);
torque = zeros(cnt, n);
t_d = zeros(cnt, n);
pos_error = zeros(cnt, 1);
rot_error = zeros(cnt, 1);
phi = zeros(cnt,1);
for i = 1 : cnt
    Td = p(t(i));
    T = forward_kin_general(robot, y(i, 1:n)) ;
    phi(i) = cal_kuka_kesai(y(i, 1:n));
    pos_error(i) = norm(Td(1:3,4) - T(1:3,4));
    rot_error(i) = norm(logR(T(1:3,1:3)' * Td(1:3,1:3)));
    torque(i,:) = controller(t(i), y(i,:)', zeros(6,n));
    t_d(i,:) = y(i,2*n +1:3*n) + y(i,n+1:2*n) * Y';
end
disp(pos_error(end));

fs = 20;
ls = 14;
tl = tiledlayout(2,3);
nexttile
plot(t,y(:,n+1:2*n),'LineWidth',2);
grid on;
ylabel('$\dot{q}$/(rad/s)', 'interpreter','latex');
xlabel("$t$/s", 'interpreter','latex', 'FontSize', fs);
set(gca,'FontSize', fs);
lg = legend('关节1', '关节2','关节3','关节4','关节5','关节6','关节7', 'Orientation','horizontal');
fontsize(lg,ls,'points')

nexttile
plot(t,t_d,'LineWidth',2);
grid on;
ylabel('$\hat{\tau}_d$/(Nm)', 'interpreter','latex');
xlabel("$t$/s", 'interpreter','latex', 'FontSize', fs);
set(gca,'FontSize', fs);
lg = legend('关节1', '关节2','关节3','关节4','关节5','关节6','关节7', 'Orientation','horizontal');
fontsize(lg,ls,'points')

nexttile
plot(t,torque,'LineWidth',2);
grid on;
ylabel('$\tau$/(Nm)', 'interpreter','latex');
xlabel("$t$/s", 'interpreter','latex', 'FontSize', fs);
set(gca,'FontSize', fs);
lg = legend('关节1', '关节2','关节3','关节4','关节5','关节6','关节7', 'Orientation','horizontal');
fontsize(lg,ls,'points')


nexttile
plot(t,pos_error*1e3,'LineWidth',2);
grid on;
ylabel('$||p_e||$/(mm)', 'interpreter','latex');
xlabel("$t$/s", 'interpreter','latex', 'FontSize', fs);
set(gca,'FontSize', fs);


nexttile
plot(t,rot_error,'LineWidth',2);
grid on;
ylabel('$||r_e||$/(rad)', 'interpreter','latex');
xlabel("$t$/s", 'interpreter','latex', 'FontSize', fs);
set(gca,'FontSize', fs);

nexttile
plot(t,phi,'LineWidth',2);
grid on;
ylabel('arm angle/(rad)', 'interpreter','latex');
xlabel("$t$/s", 'interpreter','latex', 'FontSize', fs);
set(gca,'FontSize', fs);

tl.TileSpacing = 'tight';
tl.Padding = 'compact';
% title(tl,'仿真结果', 'FontSize', fs);
set(gcf, 'Position', get(0, 'Screensize'));
% xlabel(tl, "$t$/s", 'interpreter','latex', 'FontSize', fs);
% set(gcf,'Position',[0 0 1800 1000]);
end

function [tau, daux] = PD_controller(robot, desired_motion, Kx, Bx, Kn, Bn, kesai, DOB, t, y)
n = robot.dof;
q = y(1:n);
dq = y(n + 1 : 2 * n);
[M, C, G, Jb, dJb, dM, dT, T] = m_c_g_matrix(robot,q,dq);
R = T(1:3,1:3);
p = T(1:3,4);
Vb = Jb * dq;
wb = Vb(1:3);
v = R * Vb(4:6);
dx = [wb; v];
Jh = [eye(3), zeros(3); zeros(3), R];
dJh = [zeros(3), zeros(3); zeros(3), dT(1:3,1:3)];
dJ = dJh * Jb + Jh * dJb;
J = Jh * Jb;

[Td, vel, acc] = desired_motion(t);
Rd = Td(1:3, 1:3);
pd = Td(1:3, 4);
wd = R' * vel(1:3);
vd = vel(4:6);
alphad = R' * acc(1:3);
ad = acc(4:6);
dxd = [wd;vd];
ddxd = [alphad - cross(wb, wd) ;ad];
xe = [logR(R'*Rd)'; pd - p];
dxe = dxd - dx;

% ddxc = ddxd + A_x_inv(J, M) * ((Mu_x(J, M, dJ, C) + Bx) * dxe + Kx * xe); % PD+
% if choice == 1 % pd
    ddxc = ddxd + Kx * xe + Bx * dxe;
% elseif choice == 2 % pd+
    % ddxc = ddxd + A_x_inv(J, M) * ((Mu_x(J, M, dJ, C) + Bx) * dxe + Kx * xe); % PD+
% else % passivity
    % s = dxe + Kx  * xe;
    % ddxc = ddxd + Kx * dxe + A_x_inv(J, M) * ((Mu_x(J, M, dJ, C) + Bx)) * s;
% end

a1 = pinv_J_x(J, M, ddxc - dJ * dq);

q0 = inverse_kin_kuka_robot_kesai_near(robot, Td, kesai, q);
qe = q0 - q;
qed = -dq;
a2 = null_proj(J, M, M \ (Bn * qed + Kn * qe));
% Z = null_z(J);
% dZ = derivative_null_z(J, dJ);
% ax2 = A_v(Z, M) \ ((Mu_v(Z, M, dZ, dM, C) + Bn(1)) * (-pinv_Z(Z, M) * dq) + Z' * Kn(1) * qe);
% a2 = Z * (ax2 - derivative_pinv_Z(Z, M, dZ, dM) * dq);

tau = M * (a1 + a2) + C * dq + G;

[tau, daux] = DOB(robot, tau, M, C, G, J, y);
end

function F = Wrench(t, y, robot)
F = zeros(6, robot.dof);
if t > 1 
    F(:,4) = [0, 0, 0, 0, 0, 10]';
    F(:,7) = [0, 0, 0, 0, 0, 10]';
end
end

