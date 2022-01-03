function [Mq, C, dJ] = mass_and_c_matrix2(robot, q, qd)
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
n = robot.dof;
J = zeros(6, n, n);
Mq = zeros(n, n);
dJ = zeros(6, n, n,n);
dMq = zeros(n, n, n);
C = zeros(n, n);
G = zeros([6, 6, n]);
for i = n : -1 : 1
   G(:,:,i) = [inertia(:,:,i), zeros(3);zeros(3), mass(i) * eye(3)];
   T = eye(4);
   dT = zeros([4,4,n]);
   for j = i : -1: 1
       for k = 1 : i
           [dAdT, AdT] = derivative_adjoint_T(T, dT(:,:,k));
           dJ(:,j,i,k)= dAdT * A(j,:)';
       end
       J(:,j,i) = AdT * A(j,:)';
       tform = exp_twist(-A(j,:) * q(j)) * tform_inv(M(:,:,j));
       for k = 1 : i
           if k == j
               dT(:,:,k) = (dT(:,:,k)  + T * se_twist(-A(j,:))) * tform;
           else
               dT(:,:,k) = dT(:,:,k) * tform;
           end
       end
       T = T * tform;
   end
   Mq = Mq + J(:,:,i)'*G(:,:,i)*J(:,:,i);
end
for i = 1 : n
    for j = 1 : n
        dMq(:,:,i) = dMq(:,:,i) + dJ(:,:,j,i)'*G(:,:,j)*J(:,:,j) + J(:,:,j)'*G(:,:,j)*dJ(:,:,j,i);
    end
end
for k = 1 : n
    for j = 1 : n
        for i = 1 : n
            C(k, j) = C(k, j) + 0.5 * (dMq(k,j,i) + dMq(k,i,j) - dMq(i,j,k))*qd(i);
        end
    end
end