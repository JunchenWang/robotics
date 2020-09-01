function J = derivative_q_r(r)
r = r(:);
theta = norm(r);
if theta == 0 %% deal with theta ---> 0
    J = [0, 0, 0; 0.5 * eye(3)];
    return;
end
k = sin(theta / 2) / theta;
lambda = sin(theta / 2) / theta^3 - cos(theta / 2) / (2 * theta^2);
v = k * r;
J = [-v(:)' / 2; k * eye(3) - lambda * (r * r')];

% theta_r = derivative_Nx(r);
% Jk_r = derivative_unitvector_x(r);
% k = r / norm(r);
% J2 = [-0.5 * sin(theta / 2) * theta_r
%     0.5 * cos(theta / 2) * theta_r * k(1) + sin(theta / 2) * Jk_r(1,:)
%     0.5 * cos(theta / 2) * theta_r * k(2) + sin(theta / 2) * Jk_r(2,:)
%     0.5 * cos(theta / 2) * theta_r * k(3) + sin(theta / 2) * Jk_r(3,:)];
% J2 - J
    

