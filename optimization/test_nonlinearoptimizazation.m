function test_nonlinearoptimizazation
x0 = [1, -0.1, -0.1]';
x = LevenbergMarquardt(@FletcherPowell, x0);
opt = optimset(@lsqnonlin);
opt = optimset(opt, 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'Jacobian', 'on', 'TolFun', 1e-12, 'TolX', 1e-12);%
[P,err] = lsqnonlin(@FletcherPowell, x0, [], [], opt);



function [f J] = obj(x)
f = x;
J = eye(length(x));