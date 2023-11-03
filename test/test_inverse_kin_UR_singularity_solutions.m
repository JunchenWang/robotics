function  test_inverse_kin_UR_singularity_solutions
param = [0.425, 0.3922, 0.1625, 0.1333, 0.0997, 0.0996];
port = udpport("byte");
q = [0   -1.0472    1.5359   -0.4538    0         0];
Td = forward_kin_UR(q, param);
rate = rateControl(5);
tol = [1e-5, 1e-5];
for offset = linspace(0, pi, 6)
    [qd, flag] = inverse_kin_UR(Td, param, [1, -1, -1], tol, offset);
    if flag
        set_joints(port, qd);
    else
        disp('no solution');
    end
    waitfor(rate);
end