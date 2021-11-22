function yd = computed_torque_control_law(dynamics, robot, q, qd, qdd, t, y)
persistent sum;
persistent pre_t;
if isempty(sum)
    sum = 0;
    pre_t = 0;
end
n = robot.dof;
kp = 2;
kd = 4;
ki = 0;
Kp = kp * eye(n);
Kd = kd * eye(n);
Ki = ki * eye(n);
qe = q(t) - y(1:n);
qde = qd(t) - y(n+1:end);
sum = sum + (t - pre_t) * qe;
pre_t = t;
b = qdd(t) + Kp * qe + Ki * sum +  Kd * qde;
tao = mass_matrix(robot, y(1:n)) * b + gravity_velocity_torque(robot, y(1:n), y(n+1:end));
yd = dynamics(t, y, @(t, y) tao);




