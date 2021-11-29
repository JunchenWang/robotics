function tao = computed_torque_control(robot, qd, k, t, y)
% qd is desired motion (pos, vel, acc)
persistent sum;
persistent pre_t;
if isempty(sum)
    sum = 0;
    pre_t = 0;
end
n = robot.dof;
kp = k(1);
ki = k(2);
kd = k(3);
Kp = kp * eye(n);
Kd = kd * eye(n);
Ki = ki * eye(n);
t_qd = qd(t);
qe = t_qd(1:n) - y(1:n);
qde = t_qd(n+1:2 * n) - y(n+1:end);
sum = sum + (t - pre_t) * qe;
pre_t = t;
b = t_qd(2*n+1:end) + Kp * qe + Ki * sum +  Kd * qde;
tao = mass_matrix(robot, y(1:n)) * b + gravity_velocity_torque(robot, y(1:n), y(n+1:end));




