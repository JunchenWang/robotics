function J = derivative_r_q(q)
% q(1) must be less equal than 1
q = q(:);
a = q(1);
if a == 1 %% deal with a ---> 1
    J = [-q(2:4), 2 * eye(3)];
    return;
end
if a >= 0
    v = q(2 : 4);
    t = acos(a);
    st = sqrt(1 - a^2); % sin(t)
    J = [(2 * t * a - 2 * st) / st^3 * v, 2 * t / st * eye(3)];
else
    J = -derivative_r_q(-q);
end