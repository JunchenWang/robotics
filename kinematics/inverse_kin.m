function config = inverse_kin(M, B, config, T, tol)
% M: 4x4, ef frame with respect to base
% B: nx6, screw axis of joint 1 2 ... n in body
% config: n-joint angle vector
twist = logT(tform_inv(forward_kin(M, B, config, 'b')) * T);
while norm(twist(1:3)) > tol(1) || norm(twist(4:6)) > tol(2)
    J = jacobian(B, config, 'b');
    delta = pinv(J) * twist';
    config = config + delta;
    twist = logT(tform_inv(forward_kin(M, B, config, 'b')) * T);
end

