function  test_inverse_kin_UR_solutions
param = [0.425, 0.3922, 0.1625, 0.1333, 0.0997, 0.0996];
port = udpport("byte");
% q = query_joints(port);
% q = [0   -1.0472    1.5359   -0.4538    1.5708         0];
q = [0   -1.0472    1.5359   -0.4538    0         0];
Td = forward_kin_UR(q, param);
disp(Td);
% Td(1:3,1:3) = Td(1:3,1:3) * rotz(45);
tol = [1e-5,1e-5];
rate = rateControl(1);
for cfg1 = -1 : 2 : 1
    for cfg2 =  -1 : 2 : 1
        for cfg3 = -1 : 2 : 1
            cfg = [cfg1, cfg2, cfg3];
            [qd, flag] = inverse_kin_UR(Td, param, cfg, tol);
            if flag == 1 
                disp(qd');
                set_joints(port, qd);
                waitfor(rate);
            else
                disp(cfg);
                disp('no solution');
            end
        end
    end
end