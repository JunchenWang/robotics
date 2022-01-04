function tao = computed_torque_control_task(robot, Td, w_v_d, alpha_a_d, k, freq, t, y)
% Xd is desired motion in task space
persistent sum;
persistent pre_t;
persistent pre_Xe;
if isempty(sum)
    sum = zeros(6,1);
    pre_t = t;
    pre_Xe = zeros(6,1);
end
cnt = round(t * freq + 1);
n = robot.dof;
q = y(1:n);% + 1e-3 * rand(1, 7); % noise
qd = y(n + 1:end);% + 1e-3 * rand(1, 7); %noise
kp = k(1);
ki = k(2);
kd = k(3);
factor = k(4);
Kp = kp * [factor * eye(3), zeros(3); zeros(3), eye(3)];
Kd = kd * [factor * eye(3), zeros(3); zeros(3), eye(3)];
Ki = ki * [factor * eye(3), zeros(3); zeros(3), eye(3)];
Xd = Td(:,:,cnt);
Rd = Xd(1:3, 1:3);
wd = Rd' * w_v_d(1:3, cnt);
vd = Rd' * w_v_d(4:6, cnt);
alphad = Rd' * alpha_a_d(1:3, cnt);
ad = Rd' * alpha_a_d(4:6, cnt);
V_d = [wd;vd];
dXd = Xd * se_twist(V_d);
dV_d = [alphad;ad - cross(wd, vd)];
[dJb, Jb, dT, T] = derivative_jacobian_matrix(robot, q, qd);
[dInvT, invT] = derivative_tform_inv(T, dT);
Dist = invT * Xd;
Xe = logT(Dist)';
dDist = dInvT * Xd + invT * dXd;
[dAdT, AdT] = derivative_adjoint_T(Dist, dDist);
se_V = invT * dT;
V = zeros(6,1);
V(1:3) = skew_mat2vec(se_V(1:3,1:3));
V(4:6) = se_V(1:3, 4);
Ve = AdT * V_d - V;
sum = sum + (t - pre_t) * (Xe + pre_Xe) / 2;
a = dAdT * V_d + AdT * dV_d + Kp*Xe + Ki*sum + Kd * Ve;
Mq = mass_matrix(robot, y(1:n));
ret = lsqminnorm(Jb, [a, V, dJb]);%一次性求解
qd = ret(:, 2);
Aqa = Mq * ret(:, 1);
MqJb_invdJb = Mq * ret(:, 3:end);
% qd = lsqminnorm(Jb, V); %%非常重要！！！！不能使用当前状态的qd
% Jb_invdJb = lsqminnorm(Jb, dJb);
% Jb_inv = pinv(Jb);
% Aq =  Mq * Jb_inv;
hqqd = gravity_velocity_torque(robot, q, qd);
etaqv =  hqqd - MqJb_invdJb * qd;
tao = Aqa + etaqv;
pre_t = t;
pre_Xe = Xe;
% disp(tao);




