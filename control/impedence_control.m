function [re, pe] = impedence_control(robot, Xd, Vd, M, B, K, y, f, dt)
n = robot.dof;
q = y(1:n);
qd = y(n + 1 : 2 * n);
% velocity noise
% qd = max(qd) * randn(n,1); 
[Jb, X] = jacobian_matrix(robot, q);
R = X(1:3,1:3);
p = X(1:3,4);
Rd = Xd(1:3, 1:3);
pd = Xd(1:3,4);
V = Jb * qd;
pe = R' * (pd - p);
ped = -so_w(V(1:3)) * pe + R' * Vd(4:6) - V(4:6);
re = logR(R'*Rd)';

norm_r = norm(re);
if norm_r < 1e-12
    Ar = eye(3);
else
    sr = vec2skew_mat(re);
    Ar = eye(3) - (1 - cos(norm_r)) / norm_r^2 * sr + (norm_r - sin(norm_r))/norm_r^3 * sr*sr;
end
% disp(norm(Ar - eye(3)));
red = Ar \ (Rd' * Vd(1:3) - Rd'*R*V(1:3));
% red =  V(1:3) - R'* Vd(1:3);
pedd = M(:,:,2) \ (f(4:6) - B(:,:,2) * ped - K(:,:,2) * pe);
redd = M(:,:,1) \ (f(1:3) - B(:,:,1) * red - K(:,:,1) * re);

% pe = pe + ped * dt + 0.5 * pedd * dt^2;
ped = ped + dt * pedd;
pe = pe + dt * ped;

% re = re + red * dt + 0.5 * redd * dt^2;
red = red + dt * redd;
re = re + dt * red;

% R = Rd * exp_w(re)';
% p = pd - R * pe;
% 
% X = [R, p; 0, 0, 0, 1];

