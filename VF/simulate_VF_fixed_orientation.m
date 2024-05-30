function simulate_VF_fixed_orientation

port = udpport("byte");
robot = convert_robot_tree2(importrobot('urdf\iiwa7\iiwa7.urdf'));
nominal_robot = robot;
% nominal_robot.mass = 1.2*robot.mass;% error
n = robot.dof;
Kx = 3000 * eye(5);
Bx = 300 * eye(5);
Bn = 2 * eye(n);
tspan = [0, 5];
MassMatrix = @(t, y) [eye(n), zeros(n, 2 * n); zeros(n), mass_matrix(robot, y(1:n)), zeros(n); zeros(n, 2 * n), eye(n)];
opts = odeset('Mass',MassMatrix,'OutputFcn',@(t, y, flag) odeplot_micsys(t, y, flag, port, robot));

y0 = zeros(3*n,1);
y0(1:n) = [0 75 0 -94 0 -81 0] / 180 * pi;
ptp_move(port, y0(1:n)');
Ts = forward_kin_general(robot, y0);
r0 = logR(Ts(1:3,1:3));
x0 = [r0(:); Ts(1:2,4)];
p = @(t) desired_task(t, x0);

Y = 100 * eye(n);
DOB = @(robot, tau, M, C, G, J, y) manipulator_DOB(robot, Y, tau, M, C, G, J, y);
controller = @(t, y, Fext) VF_controller(nominal_robot, p, Kx, Bx, Bn, DOB, t, y);
target_sysm = @(t, y) manipulator_dynamics_general(robot, controller, @(t, y) Wrench(t, y, robot), t, y);
[t,y] = ode15s(target_sysm,tspan,y0,opts);
cnt = length(t);
torque = zeros(cnt, n);
t_d = zeros(cnt, n);
rot_error = zeros(cnt, 3);
x_error = zeros(cnt, 2);
phi = zeros(cnt,1);
for i = 1 : cnt
    xd = p(t(i));
    T = forward_kin_general(robot, y(i, 1:n)) ;
    x = logR(T(1:3, 1:3))';
    phi(i) = cal_kuka_kesai(y(i, 1:n));
    rot_error(i,:) = x - xd(1:3);
    x_error(i,:) = T(1:2,4) - xd(4:5);
    torque(i,:) = controller(t(i), y(i,:)', zeros(6,n));
    t_d(i,:) = y(i,2*n +1:3*n) + y(i,n+1:2*n) * Y';
end
% disp(pos_error(end));

fs = 20;
ls = 14;
tl = tiledlayout(2,3);
nexttile
plot(t,y(:,n+1:2*n),'LineWidth',2);
grid on;
ylabel('$\dot{q}$/(rad/s)', 'interpreter','latex');
set(gca,'FontSize', fs);
lg = legend('J1', 'J2','J3','J4','J5','J6','J7', 'Orientation','horizontal');
fontsize(lg,ls,'points')

nexttile
plot(t,t_d,'LineWidth',2);
grid on;
ylabel('$\tau_d$/(Nm)', 'interpreter','latex');
set(gca,'FontSize', fs);
lg = legend('J1', 'J2','J3','J4','J5','J6','J7', 'Orientation','horizontal');
fontsize(lg,ls,'points')

nexttile
plot(t,torque,'LineWidth',2);
grid on;
ylabel('$\tau$/(Nm)', 'interpreter','latex');
set(gca,'FontSize', fs);
lg = legend('J1', 'J2','J3','J4','J5','J6','J7', 'Orientation','horizontal');
fontsize(lg,ls,'points')


nexttile
plot(t,rot_error,'LineWidth',2);
grid on;
ylabel('$x_e$/(rad)', 'interpreter','latex');
set(gca,'FontSize', fs);


nexttile
plot(t,x_error * 1e3,'LineWidth',2);
grid on;
ylabel('$x$/(mm)', 'interpreter','latex');
set(gca,'FontSize', fs);

nexttile
plot(t,phi,'LineWidth',2);
grid on;
ylabel('arm angle/(rad)', 'interpreter','latex');
set(gca,'FontSize', fs);

tl.TileSpacing = 'tight';
tl.Padding = 'compact';
title(tl,'Simulation Result', 'FontSize', fs);
xlabel(tl,"$t$/s", 'interpreter','latex', 'FontSize', fs);
set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf,'Position',[0 0 1800 1000]);
end

function [xd, dxd, ddxd] = desired_task(t, x0)
    
     xd = x0(:);
     dxd = zeros(5,1);
     ddxd = zeros(5,1);

end

function [tau, daux] = VF_controller(robot, desired_motion, Kx, Bx, Bn, DOB, t, y)
n = robot.dof;
q = y(1:n);
dq = y(n + 1 : 2 * n);
[M, C, G, Jb, dJb, dM, dT, T] = m_c_g_matrix(robot,q,dq);

[J, dJ, x, dx] = orientation_jacobian(Jb, dJb, T, dT);


[xd, dxd, ddxd] = desired_motion(t);

xe = xd - x;
dxe = dxd - dx;

ddxc = ddxd + A_x_inv(J, M) * ((Mu_x(J, M, dJ, C) + Bx) * dxe + Kx * xe); % PD+
% if choice == 1 % pd
%     ddxc = ddxd + Kx * xe + Bx * dxe;
% elseif choice == 2 % pd+
%     ddxc = ddxd + A_x_inv(J, M) * ((Mu_x(J, M, dJ, C) + Bx) * dxe + Kx * xe); % PD+
% else % passivity
%     s = dxe + Kx  * xe;
%     ddxc = ddxd + Kx * dxe + A_x_inv(J, M) * ((Mu_x(J, M, dJ, C) + Bx)) * s );
% end

a1 = pinv_J_x(J, M, ddxc - dJ * dq);

qed = -dq;
a2 = null_proj(J, M, M \ (Bn * qed));
% Z = null_z(J);
% dZ = derivative_null_z(J, dJ);
% ax2 = A_v(Z, M) \ ((Mu_v(Z, M, dZ, dM, C) + Bn(1)) * (-pinv_Z(Z, M) * dq) + Z' * Kn(1) * qe);
% a2 = Z * (ax2 - derivative_pinv_Z(Z, M, dZ, dM) * dq);

tau = M * (a1 + a2) + C * dq + G;

[tau, daux] = DOB(robot, tau, M, C, G, J, y);
end

function F = Wrench(t, y, robot)
F = zeros(6, robot.dof);
if t > 2 && t < 3
    F(:,7) = [0, 1, 0, -5, 0, 0]';
    F(:,4) = [0, 0, 0, 0, 0, 5]';
    % F(:,7) = [0, 0, 0, 0, 0, 10]';
elseif t > 3 && t < 4
    F(:,7) = [0, -1, 0, 5, 0, 0]';
    F(:,4) = [0, 0, 0, 0, 0, -5]';
end
end

