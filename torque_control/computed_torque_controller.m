function tao = computed_torque_controller(robot, Td, w_v_d, alpha_a_d, k, freq, t, y)
% Xd is desired motion in task space
persistent sum;
persistent pre_t;
persistent pre_xe;
if isempty(sum)
    sum = zeros(6,1);
    pre_t = t;
    pre_xe = zeros(6,1);
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
if size(Td, 3) > 1
    Xd = Td(:,:,cnt);
    Rd = Xd(1:3, 1:3);
    %% 期望姿态，速度和加速度转化到body frame下
    wd = Rd' * w_v_d(1:3, cnt);
    vd = Rd' * w_v_d(4:6, cnt);
    alphad = Rd' * alpha_a_d(1:3, cnt);
    ad = Rd' * alpha_a_d(4:6, cnt);
else
    Xd = Td;
    Rd = Xd(1:3, 1:3);
    %% 期望姿态，速度和加速度转化到body frame下
    wd = Rd' * w_v_d(1:3);
    vd = Rd' * w_v_d(4:6);
    alphad = Rd' * alpha_a_d(1:3);
    ad = Rd' * alpha_a_d(4:6);
end
Vd = [wd;vd];
dVd = [alphad;ad - cross(wd, vd)];
dXd = Xd * se_twist(Vd);
%% 实际姿态，速度在bodya frame下
[dJb, Jb, dX, X] = derivative_jacobian_matrix(robot, q, qd);
[dInvX, invX] = derivative_tform_inv(X, dX);
Xe = invX * Xd;
xe = logT(Xe)';
dXe = dInvX * Xd + invX * dXd;
[dAdXe, AdXe] = derivative_adjoint_T(Xe, dXe);
V = Jb * qd;
Ve = AdXe * Vd - V;
sum = sum + (t - pre_t) * (xe + pre_xe) / 2;
a = dAdXe * Vd + AdXe * dVd + Kp * xe + Ki*sum + Kd * Ve;
a = lsqminnorm(Jb, a - dJb * qd);
tao = inverse_dynamics(robot, q, qd, a, zeros(6,1));
pre_t = t;
pre_xe = xe;
% disp(tao);




