function simulate_robot_dynamics
%% read robot parameters
port = udpport("byte");
robot = readRobotJson('iiwa7.json');
n = robot.dof;
MassMatrix = @(t, y) [eye(n), zeros(n); zeros(n), mass_matrix(robot, y(1:n))];
opts = odeset('Mass',MassMatrix,'OutputFcn',@(t, y, flag) odeplot(t, y, flag, port, robot));
%% initial state
y0 = zeros(2*n,1);
y0(1:n) = 1e-5*rand(n,1);
%% ptp move
ptp_move(port, y0(1:n));

%% simulate
tau = @(t, y) controller(robot, t, y);
fext = @(t, y) Wrench(t, y, robot);
target = @(t, y) manipulator_dynamics(robot, tau, fext, t, y);
tspan = [0, 2];
[t,y] = ode45(target,tspan,y0,opts);
q = y(:,1:n); % joint pos
dq = y(:, n+1:2*n); % joint velocity
%% plot
fs = 24;
ls = 24;
tl = tiledlayout(2,1);
nexttile
plot(t,q, 'LineWidth',2);
ylabel('$q$/ rad', 'interpreter','latex');
set(gca,'FontSize', fs);
lg = legend('J1', 'J2','J3','J4','J5','J6','J7', 'Orientation','horizontal');
fontsize(lg,ls,'points')
% xticks([0,1,2,3,4,5,6,7,8,9,10]);

nexttile
plot(t,dq,'LineWidth',2);
grid on;
ylabel('$\dot{q}$/(rad/s)', 'interpreter','latex');
set(gca,'FontSize', fs);
lg = legend('J1', 'J2','J3','J4','J5','J6','J7', 'Orientation','horizontal');
fontsize(lg,ls,'points')

tl.TileSpacing = 'tight';
tl.Padding = 'compact';
title(tl,'Simulation Result', 'FontSize', fs);
xlabel(tl,"$t$/s", 'interpreter','latex', 'FontSize', fs);

end

%% controller
function tau = controller(robot, t, y)
    n = robot.dof;
    tau = zeros(n, 1);
end

%% dynamics
function yd = manipulator_dynamics(robot, tau, fext, t, y)
% fext is a force function applied to the robot
% fext(t,y)(:,end) is applied to TCP and others are applied to link frames
n = robot.dof;
yd = zeros(2*n,1);
q = y(1:n);
dq = y(n + 1 : 2 * n);
hqqd = gravity_velocity_torque(robot, q, dq);
yd(1:n) = dq;
ext_torque = get_ext_torque(robot, q, fext(t, y));
yd(n + 1 : 2 * n) = tau(t,y) - hqqd + ext_torque;
end

%% external force
function F = Wrench(t, y, robot)
F = zeros(6, robot.dof);
end

%% send data to visualize
function ret = odeplot(t, y, flag, port, robot)
if strcmp(flag, 'init') == 1
elseif isempty(flag)
    set_joints(port, y(1:robot.dof, end));
else
end
ret = 0;
end
