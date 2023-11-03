function  test_inverse_kin_kuka_solutions
port = udpport("byte");
% q = query_joints(port);
q = [0    0.8154   -0.1637   -1.1608         0    0.8032   -1.1713]';
set_joints(port, q);
% q = [0   -1.0472    1.5359   -0.4538    1.5708         0];
% q = [0   -1.0472    1.5359   -0.4538    0         0];
Td = forward_kin_kuka(q);
disp(Td);
% Td(1:3,1:3) = Td(1:3,1:3) * rotz(45);
tol = [1e-5,1e-5];
rate = rateControl(15);
% for cfg1 = -1 : 2 : 1
%     for cfg2 =  -1 : 2 : 1
%         for cfg3 = -1 : 2 : 1
            cfg = [1, -1, 1];
            ref = q;
            for kesai = linspace(-pi/3, pi/3, 100)
            [qd, flag] = inverse_kin_kuka_kesai(Td, cfg, kesai, tol, ref);
            ref = qd;
            if flag == 1 
                disp(qd');
                set_joints(port, qd);
                waitfor(rate);
            else
                disp(cfg);
                disp('no solution');
            end
            end
%         end
%     end
% end