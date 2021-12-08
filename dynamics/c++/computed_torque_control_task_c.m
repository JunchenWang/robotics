function tao = computed_torque_control_task_c(robot, Td, w_v_d, alpha_a_d, kp, kd, ki, t, y)
% Xd is desired motion in task space
persistent sum;
persistent pre_t;
persistent pre_Xe;
if isempty(sum)
    sum = zeros(6,1);
    pre_t = t;
    pre_Xe = zeros(6,1);
end
n = robot.dof;
q = y(1:n);
qd = y(n + 1:end);
Kp = diag(kp);
Kd = diag(kd);
Ki = diag(ki);
Xd = Td;
Rd = Xd(1:3, 1:3);
wd = Rd' * w_v_d(1:3);
vd = Rd' * w_v_d(4:6);
alphad = Rd' * alpha_a_d(1:3);
ad = Rd' * alpha_a_d(4:6);
V_d = [wd;vd];
dXd = Xd * se_twist(V_d);
dV_d = [alphad;ad - cross(wd, vd)];
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
Ve = AdT * V_d - V;
sum = sum + (t - pre_t) * (Xe + pre_Xe) / 2;
a = dAdT * V_d + AdT * dV_d + Kp*Xe + Ki*sum + Kd * Ve;
Mq = mass_matrix(robot, y(1:n));
Jb_inv = pinv(Jb);
qd = Jb_inv * V;
Aq = Mq * Jb_inv;
hqqd = velocity_torque(robot, q, qd);% 不补偿重力
etaqv =  hqqd - Aq * dJb * qd;
tao = Aq * a + etaqv;
pre_t = t;
pre_Xe = Xe;





