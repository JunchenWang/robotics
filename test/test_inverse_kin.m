function angles = test_inverse_kin
robot = UR5e;
N = 1000;
angles = zeros(N, 6);
cnt = 1;
tic;
for i = 1 : N
    q = rand(1,6);
%     q(3) = deg2rad(10) + rand * deg2rad(120);
     q(5) = deg2rad(5) + rand * deg2rad(140);
    [ang, flag] = inverse_kin_general_mex(robot, forward_kin_general(robot, q), q + (rand * 2 - 1) * 0.1, [1e-5, 1e-5]);
    if flag == 0
        angles(cnt,:) = q;
        cnt = cnt + 1;
        disp(q);
    end
end
angles(cnt:end,:) = [];
toc;