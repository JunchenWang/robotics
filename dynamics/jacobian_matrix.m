function [Jb, T] = jacobian_matrix(A, M, ME, q)
T = ME;
for i = 7 : -1 : 1
    Jb(:,i) = adjoint_T(tform_inv(T)) * A(i,:)';
    T = M(:,:,i) * exp_twist(A(i,:) * q(i)) * T;
end