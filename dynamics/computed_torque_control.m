function tao = computed_torque_control(robot, q, qd, qdd, kp, ki, kd, freq, t, y)
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
cnt = round(freq * t + 1);
Kp = diag(kp);
Kd = diag(kd);
Ki = diag(ki);
qe = q(:,cnt) - y(1:n);
qde = qd(:, cnt) - y(n+1:end);
sum = sum + 0.5 * (pree + qe) * (t - pre_t);
pre_t = t;
pree = qe;
b = qdd(:,cnt) + Kp * qe + Ki * sum +  Kd * qde;
tao = mass_matrix(robot, y(1:n)) * b + gravity_velocity_torque(robot, y(1:n), y(n+1:end));




