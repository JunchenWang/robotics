function tao = computed_torque_control_task(robot, Td, w_v_d, alpha_a_d, k, freq, t, y)
% Xd is desired motion in task space
persistent sum;
persistent pre_t;
if isempty(sum)
    sum = 0;
    pre_t = 0;
end

cnt = round(t * freq + 1);
n = robot.dof;
q = y(1:n);
qd = y(n + 1:end);
kp = k(1);
ki = k(2);
kd = k(3);
Kp = kp * eye(6);
Kd = kd * eye(6);
Ki = ki * eye(6);
Xd = Td(:,:,cnt);
Rd = Xd(1:3, 1:3);
wd = Rd' * w_v_d(1:3, cnt);
vd = Rd' * w_v_d(4:6, cnt);
alphad = Rd' * alpha_a_d(1:3, cnt);
ad = Rd' * alpha_a_d(4:6, cnt);
twist_d = [wd;vd];
dXd = Xd * se_twist(twist_d);
dtwist_d = [alphad;ad - cross(wd, vd)];
[Jb, T, dJb, dT] = derivative_jacobian_matrix(robot, q, qd);
[dInvT, invT] = derivative_tform_inv(T, dT);
Dist = invT * Xd;
Xe = logT(Dist)';
dDist = dInvT * Xd + invT * dXd;
[dAdT, AdT] = derivative_adjoint_T(Dist, dDist);
se_V = invT * dT;
V = zeros(6,1);
V(1:3) = skew_mat2vec(se_V(1:3,1:3));
V(4:6) = se_V(1:3, 4);
Ve = AdT * twist_d - V;
sum = sum + (t - pre_t) * Xe;
a = dAdT * twist_d + AdT * dtwist_d + Kp*Xe + Ki*sum + Kd * Ve;
Mq = mass_matrix(robot, y(1:n));
qd = lsqminnorm(Jb, V); %%非常重要！！！！不能使用当前状态的qd
Jb_inv = pinv(Jb);
Aq =  Mq * Jb_inv;
hqqd = gravity_velocity_torque(robot, q, qd);
etaqv =  hqqd - Aq*dJb*qd;
tao = Aq * a + etaqv;
pre_t = t;





