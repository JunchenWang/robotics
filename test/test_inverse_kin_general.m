function [err1, err2, failure1, failure2, t1, t2, cnt1, cnt2] = test_inverse_kin_general(N)
if nargin < 1
    N = 10000;
end
% Rd = eul2rotm([180, 180, 180] / 180 * pi);
% pd = [-525.9610,-133.6040,654.6057]' / 1000;
robot = convert_robot_tree2(importrobot('urdf\ur_5e-calibrated\ur_description\urdf\ur5e-A302.urdf'));
% robot = convert_robot_tree2(importrobot('urdf\iiwa7\iiwa7.urdf'));
n = robot.dof;
err1 = zeros(N, 2);
err2 = zeros(N, 2);
t1 = 0;
t2 = 0;
cnt1 = 0;
cnt2 = 0;
tol = [1e-5, 1e-5];
failure1 = zeros(100, n);
failure2 = zeros(100, n);
for i =1 : N
    q = rand(n, 1) * 2 * pi - pi;
    % q = [-0.0635   -2.0865    3.0076    1.3364    0.0030   -0.1817]';
    Td = forward_kin_general(robot, q);
    % Td(1:3,1:3) = eye(3);
    ref = mod(q + (rand(n,1) * 2 - 1) / 6 + pi, 2*pi) - pi;
    % Td = [Rd, pd; 0 0 0 1];
    tic;
    [angles, flag] = inverse_kin_general(robot, Td, ref, tol);
    t = toc;
    t1 = t1 + t;
    T = forward_kin_general(robot, angles);
    err1(i, 1) = norm(logR(T(1:3,1:3)'*Td(1:3,1:3)));
    err1(i, 2) = norm(T(1:3,4)- Td(1:3,4));
    if ~flag
        cnt1 = cnt1 + 1;
        failure1(cnt1, :) = q;
    end
    % tic;
    % [angles, flag] = inverse_kin_general2(robot, Td, ref, tol);
    % t = toc;
    % t2 = t2 + t;
    % T = forward_kin_general(robot, angles);
    % err2(i, 1) = norm(logR(T(1:3,1:3)'*Td(1:3,1:3)));
    % err2(i, 2) = norm(T(1:3,4)- Td(1:3,4));
    % if ~flag
    %     cnt2 = cnt2 + 1;
    %     failure2(cnt2, :) = q;
    % end
end
% err1 = mean(err1, 1);
% err2 = mean(err2, 1);
end

