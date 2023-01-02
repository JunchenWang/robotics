function angles = inverse_kin_kuka_robot_kesai_near(robot, T, kesai, ref)
% T is in meter
angles = [];
dist = 1e9;
for cfg1 = -1 : 2 : 1
    for cfg2 =  -1 : 2 : 1
        for cfg3 = -1 : 2 : 1
            cfg = [cfg1, cfg2, cfg3];
            ang = inverse_kin_kuka_robot_kesai(robot, T, cfg, kesai, ref);
            if ~isempty(ang) && norm(ang-ref)<dist
                dist = norm(ang-ref);
                angles = ang;
            end
        end
    end
end
if isempty(angles) || limit_check_kuka(angles)
    angles = [];
end