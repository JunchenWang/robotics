function [angles, flag] =UR_inverse_kin(robot, Td, ref)
% test: q = [0.1328   -1.6864   -0.0698    0.7795    1.1255   -0.6565];
% test: q =  [-0.0635   -2.0865    3.0076    1.3364    0.0030   -0.1817];
% wrong UR_inverse_kin_cfg
dist = 1e10;
angles = zeros(1,6);
flag = 0;
for cfg1 = -1 : 2 : 1
    for cfg3 = -1 : 2 : 1
        for cfg5 = -1 : 2 : 1
            cfg = [cfg1, cfg3, cfg5];
            angles_cfg = UR_inverse_kin_cfg(Td, cfg, ref(6));
            display(angles_cfg);
            err = norm(angles_cfg - ref);
            if err < dist
                dist = err;
                angles = angles_cfg;
            end
        end
    end
end
axang = rotm2axang(Td(1:3,1:3));
rd = axang(1:3)' * axang(4);
td = Td(1:3,4);
% avoid large condition number of Ja
if abs(angles(3)) < 0.05
    angles(3) = 0.05;
end
if abs(angles(5)) < 0.05
    angles(5) = 0.05;
end
[Jb, T] = jacobian_matrix(robot, angles);
axang = rotm2axang(T(1:3,1:3));
r = axang(1:3)' * axang(4);
t = T(1:3,4);
cnt = 0;
while (norm(rd - r) > deg2rad(0.001) || norm(td - t) > 1e-4) && cnt < 6
    Ja = analytic_jacobian_matrix(Jb, T);
    delta = lsqminnorm(Ja, [rd;td] - [r;t]);
    angles = angles + delta';
    cnt = cnt + 1;
    [Jb, T] = jacobian_matrix(robot, angles);
    axang = rotm2axang(T(1:3,1:3));
    r = axang(1:3)' * axang(4);
    t = T(1:3,4);
end
if cnt < 6
    flag = 1;
end