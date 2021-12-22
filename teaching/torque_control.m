function tao = torque_control(robot, q, qd, qdd, kp, ki, kd, t, y)
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
qt = q(t);
qdt = qd(t);
qddt = qdd(t);
Kp = diag(kp);
Kd = diag(kd);
Ki = diag(ki);
qe = qt - y(1:n);
qde = qdt - y(n+1:end);
sum = sum + 0.5 * (pree + qe) * (t - pre_t);
pre_t = t;
pree = qe;
b = qddt + Kp * qe + Ki * sum +  Kd * qde;
tao = MassMatrix(robot, y(1:n)) * b + CMatrix(robot, y(1:n), y(n+1:end)) * y(n+1:end) + GMatrix(robot, y(1:n));
% tao = b;



