function [angles, flag] = inverse_kin_general2(robot, Td, ref, tol)
% test: q = [0.1328   -1.6864   -0.0698    0.7795    1.1255   -0.6565];
% test: q =  [-0.0635   -2.0865    3.0076    1.3364    0.0030   -0.1817];

[angles, flag] = inverse_kin_general_lsqnonlin(robot, Td, ref, tol);
% angles = ref;
% rd = logR(Td(1:3,1:3))';
% pd = Td(1:3,4);
% [Jb, T] = jacobian_matrix(robot, angles);
% r = logR(T(1:3,1:3))';
% p = T(1:3,4);
% cnt = 0;
% b = [r;p] - [rd;pd];
% 
% while (norm(b(1:3))  > tol(1) || norm(b(4:6)) > tol(2)) && cnt < 100
%     delta = pinv(Jb) * ([w_dr_A(r), zeros(3); zeros(3),T(1:3,1:3)'] * b);
%     angles = angles - delta';
%     cnt = cnt + 1;
%     [Jb, T] = jacobian_matrix(robot, angles);
%     r = logR(T(1:3,1:3))';
%     p = T(1:3,4);
%     b = [r;p] - [rd;pd];
% end
% angles = mod(angles + pi, 2*pi) - pi;
% if cnt < 10
%     flag = 1;
% else
%     disp('nolinear');
%     [angles, flag] = inverse_kin_general_lsqnonlin(robot, Td, ref, tol);
% end
end

function [angles, flag] = inverse_kin_general_lsqnonlin(robot, Td, ref, tol)
% tol = [1e-5, 1e-5]
fun = @(x) objFun(x, Td, robot);
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Algorithm', 'levenberg-marquardt',...
                        'FunctionTolerance',1e-10, 'Display','off','MaxIterations', 100);
[angles,~,residual, exitflag]  = lsqnonlin(fun,ref,[],[],options);
angles = mod(angles + pi, 2*pi) - pi;
if norm(residual(1:3)) < tol(1) && norm(residual(4:6)) < tol(2)
    flag = 1;
else
    flag = 0;
end
end
function [f, g] = objFun(x, Td, robot)
rd = logR(Td(1:3,1:3))';
pd = Td(1:3,4);
[Jb, T] = jacobian_matrix(robot, x);
g = analytic_jacobian_matrix(Jb, T);
r = logR(T(1:3,1:3))';
p = T(1:3,4);
f = [r;p]-[rd;pd];
end