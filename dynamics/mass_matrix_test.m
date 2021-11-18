function Mq = mass_matrix_test(robot, q)
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
n = length(mass);
J = zeros(6, n, n);
Mq = zeros(n, n);
AA = A;
for i = 1 : n
   G = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
   T = M(:,:,i) * exp_twist(A(i,:) * q(i));
   J(:, i, i)  = AA(i,:)';
   for j = 1 : i - 1
       AA(j,:) = adjoint_T(tform_inv(T)) * AA(j,:)';
       J(:, j, i)  = AA(j,:)';
   end
   Mq = Mq + J(:,:,i)'*G*J(:,:,i);
end