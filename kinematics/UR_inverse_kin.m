function [angles, flag] =UR_inverse_kin(robot, Td, ref)
dist = 1e10;
angles = zeros(1,6);
flag = 0;
for cfg1 = -1 : 2 : 1
    for cfg2 = -1 : 2 : 1
        for cfg3 = -1 : 2 : 1
            cfg = [cfg1, cfg2, cfg3];
            [angles_cfg, flag_cfg] = UR_inverse_kin_cfg(Td, cfg, ref(6));
            if flag_cfg == 0
                continue;
            end
            err = norm(angles_cfg - ref);
            if err < dist
                flag = 1;
                dist = err;
                angles = angles_cfg;
            end
        end
    end
end
axang = rotm2axang(Td(1:3,1:3));
rd = axang(1:3)' * axang(4);
td = Td(1:3,4);
[Jb, T] = jacobian_matrix(robot, angles);
axang = rotm2axang(T(1:3,1:3));
r = axang(1:3)' * axang(4);
t = T(1:3,4);
if flag == 1
    cnt = 0;
    while (norm(rd - r) > deg2rad(0.001) || norm(td - t) > 1e-4) && cnt < 5
        Ja = analytic_jacobian_matrix(Jb, T);
        delta = lsqminnorm(Ja, [rd;td] - [r;t]);
        angles = angles + delta';
        cnt = cnt + 1;
        [Jb, T] = jacobian_matrix(robot, angles);
        axang = rotm2axang(T(1:3,1:3));
        r = axang(1:3)' * axang(4);
        t = T(1:3,4);
    end
    disp(cnt);
end