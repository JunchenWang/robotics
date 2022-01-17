function [angles, flag] =UR_inverse_kin(T, ref)
dist = 1e10;
angles = zeros(1,6);
flag = 0;
delta_theta = [ -1.30565669838322851e-08, -0.588159115180872494, 0.826314252559053, -0.238165200723693349, -1.07362188606074938e-06, -2.93069160793579808e-07];
for cfg1 = -1 : 2 : 1
    for cfg2 = -1 : 2 : 1
        for cfg3 = -1 : 2 : 1
            cfg = [cfg1, cfg2, cfg3];
            [angles_cfg, flag_cfg] = UR_inverse_kin_cfg(T, cfg, ref(6) + delta_theta(6));
            if flag_cfg == 0
                continue;
            end
            angles_cfg = angles_cfg - delta_theta;
            err = norm(angles - ref);
            if err < dist
                flag = 1;
                dist = err;
                angles = angles_cfg;
            end
        end
    end
end