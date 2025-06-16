function KFTest
% DC simulator
Q = 1e-5;
R = 0.1^2;
interval = 0.05;
N = 1000;
Y = 0.5;
y = Y + sqrt(R) * randn(1,N);
t = interval * [1 : N];
A = 1;
H = 1;
x0 = 0.3;
P0 = 1;
[x P] = EKalmanFilter(@fun, @hun, x0, P0, y, R, Q);
plot(t, y, 'r');
hold on;
plot([0 t], Y * ones(1, N + 1), 'm');
% DrawUncertainty(x, P, 'c', 1, interval);
plot([0 t], x, 'g*');
hold off;

function [f Jx Jw] = fun(x, w)
f = x + w;
Jx = 1;
Jw = 1;

function [h Jx Jv] = hun(x, v)
h = x + v;
Jx = 1;
Jv = 1;