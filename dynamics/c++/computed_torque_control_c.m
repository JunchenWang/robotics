function tao = computed_torque_control_c(robot, q, qd, qdd, kp, ki, kd, t, y)
% qd is desired motion (pos, vel, acc)
n = robot.dof;
persistent sum;
persistent pre_t;
persistent pree;
if isempty(sum)
    sum = zeros(n, 1);
    pre_t = t;
    pree = zeros(n, 1);
end
Kp = diag(kp);
Kd = diag(kd);
Ki = diag(ki);
qe = q - y(1:n);
qde = qd - y(n+1:end);
sum = sum + 0.5 * (pree + qe) * (t - pre_t);
pre_t = t;
pree = qe;
b = qdd + Kp * qe + Ki * sum +  Kd * qde;
tao = mass_matrix(robot, y(1:n)) * b + velocity_torque(robot, y(1:n), y(n+1:end));