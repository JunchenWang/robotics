function BayesianEstimateTest
Q = 1e-5;
R = 0.1^2;
N = 500;
Y = 0.5;
y = Y + sqrt(R) * randn(1,N);
x0 = 0.3;
P0 = 1;
% xxx = RegularizedBayesianFilter(@BayesianHMM, @BayesianLikelihood, x0, @BayesianHMMP0, y, R, Q);
IsDebug = 0;
x = SISBayesianFilter(@BayesianHMM, @BayesianLikelihood, x0, @BayesianHMMP0, y, R, Q, IsDebug);
xx = EKalmanFilter(@fun, @hun, x0, P0, y, R, Q);

plot(0 : N, xx,'*m')
hold on;
plot(1 : N, y, 'b');
plot(0 : N, x, 'g*');
% plot(0 : N, x, 'y*');
function p = BayesianHMMP0(op, x0)
if strcmp(op, 'DS') == 1
    p = 0.3 + sqrt(1) * randn(1, x0);
else
    p = normpdf(x0, 0.3, sqrt(1));
end

function p = BayesianLikelihood(yk, xk, v)
p = normpdf(yk - xk, 0, sqrt(v));

function p = BayesianHMM(op, xk, xk_1, w)
if strcmp(op, 'DS') == 1
      p = xk_1 + sqrt(w) * randn(1, xk);
else
    p = normpdf(xk - xk_1, 0, sqrt(w));
end
function [f Jx Jw] = fun(x, w)
f = x + w;
Jx = 1;
Jw = 1;

function [h Jx Jv] = hun(x, v)
h = x + v;
Jx = 1;
Jv = 1;