function [angles, flag] = inverse_kin_general(robot, Td, ref, tol)
% tol = [2e-5, 1e-4]
fun = @(x) objFun(x, Td, robot);
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Algorithm', 'levenberg-marquardt' ,'FunctionTolerance',1e-10, 'Display','off');
[angles,~,residual]  = lsqnonlin(fun,ref,[],[],options);
angles = mod(angles + pi, 2*pi) - pi;
if norm(residual(1:3)) < tol(1) && norm(residual(4:6)) < tol(2)
    flag = 1;
else
    flag = 0;
end
end
function [f, g] = objFun(x, Td, robot)
axang = rotm2axang(Td(1:3,1:3));
rd = axang(1:3)' * axang(4);
rd2 =  -axang(1:3)' * (2*pi-axang(4));
td = Td(1:3,4);
[Jb, T] = jacobian_matrix(robot, x);
g = -analytic_jacobian_matrix(Jb, T);
axang = rotm2axang(T(1:3,1:3));
r = axang(1:3)' * axang(4);
t = T(1:3,4);
if norm(rd2 - r) > norm(rd - r)
    f = [rd;td] - [r;t];
else
    f = [rd2;td] - [r;t];
end
end