function tao = inverse_dynamics(mass, inertia, A, M, ME, q, qd, qdd, F_ME)
n = length(mass);
nu0 = zeros(6, 1);
dnu0 = [0, 0, 0, 0, 0, 9.8]';
nu = zeros(6, n);
dnu = zeros(6, n);
tao = zeros(n, 1);
for i = 1 : n
    T = exp_twist(-A(i,:)*q(i))*tform_inv(M(:,:,i));
    Map = adjoint_T(T);
    nu(:, i) = Map * nu0 + A(i,:)'*qd(i);
    nu0 = nu(:, i);
    dnu(:, i) = Map * dnu0 + adjoint_twist(nu0) * A(i,:)' * qd(i) + A(i,:)'*qdd(i);
    dnu0 = dnu(:, i);
end
F = F_ME;
T = tform_inv(ME);
for i = n : -1 : 1
    G = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
    F = adjoint_T(T)'* F + G * dnu(:,i) - adjoint_twist(nu(:,i))'*(G*nu(:,i));
    tao(i) = F'*A(i,:)';
    T =  exp_twist(-A(i,:)*q(i))*tform_inv(M(:,:,i));
end