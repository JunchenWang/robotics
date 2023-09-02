function Y = regressor(robot, q, qd, a, v)
n = robot.dof;
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
ME = robot.ME;
com = robot.com;
gravity = robot.gravity;
TCP = robot.TCP;


J = zeros(6, n, n);
ppi = zeros(10, n);
G = zeros(6, 6, n);
for i = 1 : n
    ppi(:,i) = [mass(i), mass(i) * com(i,:), inertia(1,1,i), inertia(2,2,i), inertia(3,3,i),...
                inertia(1,2,i), inertia(1,3,i), inertia(2,3,i)];
    
    G(:,:,i) = spatial_inertia_matrix(inertia(:,:,i),mass(i), com(i,:));
    T = eye(4);
   for j = i : -1: 1
       J(:,j,i) = adjoint_T(T) * A(j,:)';
       T = T * exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
   end
    V = J(:,:,i) * qd;
    disp(G(:,:,i) * V - A_matrix(V) * ppi(:,i));
end


end

function A = A_matrix(V)
w = V(1:3);
v = V(4:6);
Aw = [zeros(3,1), -so_w(v), diag(w), [w(2), w(3), 0; w(1), 0, w(3); 0, w(1), w(2)]];
Av = [v, so_w(w), zeros(3,6)];
A = [Aw; Av];
end

