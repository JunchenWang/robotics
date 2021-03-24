function J = jacobian(S, config, flag)
% S: nx6, screw axis of joint 1 2 ...
% config: n-joint angle vector
% flag: 's' for sptial jacobian, 'b' for body jacobian
n = size(S, 1);
J = zeros(6, n);
T = make_tform;
if flag == 's' % spatial jacobian
    for i = 1 : n
        J(:, i) = adjoint_T(T) * S(i, :)';
        T = T * exp_twist(S(i, :) * config(i));
    end
elseif flag == 'b' % body jacobian
    for i = n : -1 : 1
        J(:, i) = adjoint_T(T) * S(i, :)';
        T = T * exp_twist(-S(i, :) * config(i));
    end
else
    error('please specify spatial or body jacobian');
end