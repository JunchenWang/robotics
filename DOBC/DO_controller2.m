function [tao, td] = DO_controller2(robot, Td, w_v_d, alpha_a_d, Kp, Kd, Bn, Kn, kesai, freq, Y, t, y)
% Xd is desired motion in task space
cnt = round(t * freq + 1);
n = robot.dof;
q = y(1:n);% + 1e-3 * rand(1, 7); % noise
qd = y(n + 1: 2 * n);% + 1e-3 * rand(1, 7); %noise
[M, C, G, Jb, dJb, dM, dX, X] = m_c_g_matrix(robot,q,qd);
R = X(1:3,1:3);
p = X(1:3,4);
Vb = Jb * qd;
wb = Vb(1:3);
v = R * Vb(4:6);
Jh = [eye(3), zeros(3); zeros(3), R];
dJh = [zeros(3), zeros(3); zeros(3), dX(1:3,1:3)];
dJb = dJh * Jb + Jh * dJb;
Jb = Jh * Jb;
td = y(2 * n + 1 : end);
if size(Td, 3) > 1
    Xd = Td(:,:,cnt);
    Rd = Xd(1:3, 1:3);
    pd = Xd(1:3, 4);
    %% 期望姿态，速度和加速度转化到Xd下
    wd = R' * w_v_d(1:3, cnt);
    vd = w_v_d(4:6, cnt);
    alphad = R' * alpha_a_d(1:3, cnt);
    ad = alpha_a_d(4:6, cnt);
else
    Xd = Td;
    Rd = Xd(1:3, 1:3);
    pd = Xd(1:3, 4);
    %% 期望姿态，速度和加速度转化到Xd下
    wd = R' * w_v_d(1:3);
    vd = w_v_d(4:6);
    alphad = R' * alpha_a_d(1:3);
    ad = alpha_a_d(4:6);
end
q0 = inverse_kin_kuka_robot_kesai_near(robot, Xd, kesai, q)';
if isempty(q0)
    error('no inverse');
end
Vd = [wd;vd];
dVd = [alphad - cross(wb, wd) ;ad];
xe = [logR(R'*Rd)'; pd - p];
dxe = Vd - [wb;v];

Z = null_z(Jb);
dZ = derivative_null_z(Jb, dJb);


% ax1 = dVd + Kp * xe + Kd * dxe;
% ax1 = dVd + A_x_inv(Jb, M) * ((Mu_x(Jb, M, dJb, C) + Kd) * dxe + Kp * xe); % PD+
s = dxe + Kp * xe;
ax1 = dVd + Kp * dxe + A_x_inv(Jb, M) * ((Mu_x(Jb, M, dJb, C) + Kd) * s ); 
a1 = pinv_J_x(Jb, M, ax1 - dJb * qd);

qe = q0 - q;
qed = -qd;
a2 = null_proj(Jb, M, M \ (Bn * qed + Kn * qe));
% ax2 = A_v(Z, M) \ ((Mu_v(Z, M, dZ, dM, C) + Bn(1)) * (-pinv_Z(Z, M) * qd) + Z' * Kn(1) * qe);
% a2 = Z * (ax2 - derivative_pinv_Z(Z, M, dZ, dM) * qd);

% tao = M * (a1 + a2) + C * qd + G;
% td = -Y * pinv_J_x(Jb, M, s);
P = Y * qd;
tao_d = td + P;
% tao_d = tao_d - M * Z * pinv_Z(Z, M) * tao_d;
tao_d = Jb' * ((Jb * (M \ Jb'))  \ (Jb * (M \ tao_d))); % minus the component projected into the null space !
disp(tao_d);
tao = M * (a1 + a2) + C * qd + G - tao_d;
td = Y * (M \ (C * qd + G - tao - P - td));








