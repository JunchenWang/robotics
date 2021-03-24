function T = forward_kin(M, S, config, flag)
% M: 4x4, ef frame with respect to base
% S: nx6, screw axis of joint 1 2 ... n
% config: n-joint angle vector
% flag: 's' for sptial screw axis, 'b' for body screw axis
Ts = exp_twist(S .* config(:));
T = make_tform(eye(3), zeros(1,3));
for i = 1 : size(Ts, 3)
    T = T * Ts(:,:,i);
end
if flag == 's' % spacial screw axis
    T = T * M;
elseif flag == 'b' % body screw axis
    T = M * T;
else
    error('please specify spatial or body screw axis');
end
