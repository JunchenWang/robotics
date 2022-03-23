function [angles, flag] = inverse_kin_general_lsqnonlin(robot, Td, ref, tol)
% tol = [1e-5, 1e-5]
fun = @(x) objFun(x, Td, robot);
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Algorithm', 'levenberg-marquardt',...
                        'FunctionTolerance',1e-10, 'Display','off');
[angles,~,residual]  = lsqnonlin(fun,ref,[],[],options);
angles = mod(angles + pi, 2*pi) - pi;
if norm(residual(1:3)) < tol(1) && norm(residual(4:6)) < tol(2)
    flag = 1;
else
    flag = 0;
end
end
function [f, g] = objFun(x, Td, robot)
rd = logR(Td(1:3,1:3))';
norm_rd = norm(rd);
if norm_rd ~= 0
    rd2 = (2*pi - norm_rd) * (-rd / norm_rd);
else
    rd2 = rd;
end
td = Td(1:3,4);
[Jb, T] = jacobian_matrix(robot, x);
g = -analytic_jacobian_matrix(Jb, T);
r = logR(T(1:3,1:3))';
t = T(1:3,4);
if norm(rd2 - r) > norm(rd - r)
    f = [rd;td] - [r;t];
else
    f = [rd2;td] - [r;t];
end
end