function x = TrackingFiltering(y, offset)
N = size(y, 2);
xx = zeros(6, N);
for i = 1 : N
    xx(:,i) = ConcatenateFrame(y(:,i), offset);
end
x0 = xx(:,1);
cy = InverseFrame(offset);
Q = eye(6);
R = cov(y');

P0 = diag([0.001, 0.001, 0.001, 0.001, 0.001, 0.1]);
[x P] = EKalmanFilter(@fun, @hun, x0, P0, y, R, Q, [], cy);
plot(xx(6,:), 'r');
hold on;
plot(x(6,:), 'b');
hold off;

function [f Jx Jw] = fun(x, w)
f = x + w;
Jx = eye(6);
Jw = eye(6);

function [h Jx Jv] = hun(x, v, cy)
[Jx Jf1] = derivative_f2f1(x, cy);
Jv = eye(6);
h = ConcatenateFrame(x, cy) + v;
