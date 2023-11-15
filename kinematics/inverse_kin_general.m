function [angles, flag] = inverse_kin_general(robot, Td, ref, tol)
angles = ref(:);
rd = logR(Td(1:3,1:3))';
norm_rd = norm(rd);
if norm_rd ~= 0
    rd2 = (2*pi - norm_rd) * (-rd / norm_rd);
else
    rd2 = rd;
end
td = Td(1:3,4);
[Jb, T] = jacobian_matrix(robot, angles);
r = logR(T(1:3,1:3))';
t = T(1:3,4);
cnt = 0;
if norm(rd-r) < norm(rd2 - r)
    b = [rd;td] - [r;t];
else
    b = [rd2;td] - [r;t];
end
notYet = norm(b(1:3))  > tol(1) || norm(b(4:6)) > tol(2);
while notYet && cnt < 10
    %     Ja = analytic_jacobian_matrix(Jb, T);
    %     delta = lsqminnorm(Jb, [w_dr_A(r), zeros(3); zeros(3),T(1:3,1:3)'] * b);
    %     delta = Jb \ ([w_dr_A(r), zeros(3); zeros(3),T(1:3,1:3)'] * b);
    delta = pinv(Jb) * ([w_dr_A(r), zeros(3); zeros(3),T(1:3,1:3)'] * b);
    %     delta = mod(delta + pi, 2*pi) - pi;
    angles = angles + delta;
    cnt = cnt + 1;
    [Jb, T] = jacobian_matrix(robot, angles);
    r = logR(T(1:3,1:3))';
    t = T(1:3,4);
    if norm(rd-r) < norm(rd2 - r)
        b = [rd;td] - [r;t];
    else
        b = [rd2;td] - [r;t];
    end
    notYet = norm(b(1:3))  > tol(1) || norm(b(4:6)) > tol(2);
end
angles = mod(angles + pi, 2*pi) - pi;
if ~notYet
    flag = 1;
else
    % disp('nolinear');
    [angles, flag] = inverse_kin_general_lsqnonlin(robot, Td, ref(:), tol);
end
end