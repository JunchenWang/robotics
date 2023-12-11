function yd = joint_dynamics(joint, t, y, u, d)
% y: q, qd
yd = zeros(2,1);
yd(1) = y(2);
yd(2) = (u(t, y) - d(t, y) / joint.r - joint.B * y(2)) / joint.J;
end