function [angles, flag] =UR_inverse_kin(robot, Td, ref)
% deprecated, do not use
% test: q = [0.1328   -1.6864   -0.0698    0.7795    1.1255   -0.6565];
% test: q =  [-0.0635   -2.0865    3.0076    1.3364    0.0030   -0.1817];
dist = 1e10;
angles = zeros(1,6);
for cfg1 = -1 : 2 : 1
    for cfg3 = -1 : 2 : 1
        for cfg5 = -1 : 2 : 1
            cfg = [cfg1, cfg3, cfg5];
            angles_cfg = UR_inverse_kin_cfg(Td, cfg, ref(6));
            err = norm(angles_cfg - ref);
            if err < dist
                dist = err;
                angles = angles_cfg;
            end
        end
    end
end
% [angles, flag] = inverse_kin_general(robot, Td, angles, [1e-4, 2e-4]);
ref = angles;
axang = rotm2axang(Td(1:3,1:3));
rd = axang(1:3)' * axang(4);
rd2 =  -axang(1:3)' * (2*pi-axang(4));
td = Td(1:3,4);
[Jb, T] = jacobian_matrix(robot, angles);
axang = rotm2axang(T(1:3,1:3));
r = axang(1:3)' * axang(4);
t = T(1:3,4);
cnt = 0;
if norm(rd-r) < norm(rd2 - r)
    b = [rd;td] - [r;t];
else
    b = [rd2;td] - [r;t];
end
while (norm(b(1:3))  > 2e-5 || norm(b(4:6)) > 2e-4) && cnt < 10
    Ja = analytic_jacobian_matrix(Jb, T);
    delta = lsqminnorm(Ja, b);
    delta = mod(delta + pi, 2*pi) - pi;
    angles = angles + delta';
    cnt = cnt + 1;
    [Jb, T] = jacobian_matrix(robot, angles);
    axang = rotm2axang(T(1:3,1:3));
    r = axang(1:3)' * axang(4);
    t = T(1:3,4);
    if norm(rd-r) < norm(rd2 - r)
        b = [rd;td] - [r;t];
    else
        b = [rd2;td] - [r;t];
    end
end
angles = mod(angles + pi, 2*pi) - pi;
if cnt < 10
    flag = 1;
else
    [angles, flag] = inverse_kin_general(robot, Td, ref, [2e-5, 2e-4]);
end