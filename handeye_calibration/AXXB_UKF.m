function fX = AXXB_UKF(fA, fB)
fA = RVFrame2ThetaFrame(fA);
fB = RVFrame2ThetaFrame(fB);
n = size(fA, 2);
sigma_tA = diag([0.1 0.1 0.1]);
sigma_rA = diag([0.1 0.1 0.1]);
sigma_tB = diag([0.01 0.01 0.01]);
sigma_rB = diag([0.5 0.5 0.5]);
x0 = ones(6, 1);
P0 = eye(6);
Cy = cell(1, n);
Cy{1} = fB(:, 1);
y = cell(1, n);
y{1} = fA(:, 1);
R = cell(1, n);
R{1} = blkdiag(sigma_rB, sigma_tB, sigma_rA, sigma_tA);
for i = 2 : n
    Cy{i} = [Cy{i - 1}; fB(:, i)];
    y{i} = [y{i - 1}; fA(:, i)];
    R{i} = blkdiag(R{i - 1}, sigma_rB, sigma_tB, sigma_rA, sigma_tA);
end
Q = 10 * blkdiag(sigma_rB + sigma_rA, sigma_tB + sigma_tA);
[x p] = UKalmanFilter(@f, @h, x0, P0, y, R, Q, [], Cy, [0.25 3 0]);

rms = zeros(n + 1, 1);
for i = 1 : n + 1
    rms(i) = sqrt(trace(p(:, :, i)));
end
subplot(1, 2, 1);
plot(x');
subplot(1, 2, 2);
plot(rms);
fX = x(:, end);
fX = ThetaFrame2RVFrame(fX);


function yk = h(xk, vk, Cy)
n = length(Cy) / 6;
yk = zeros(6 * n, 1);
for i = 1 : n
    fxfb = ConcatenateThetaFrame(xk, Cy(6 * i - 5 : 6 * i) - vk(12 * i - 11 : 12 * i - 6));
    yk(6 * i - 5 : 6 * i) = ConcatenateThetaFrame(fxfb, InvertThetaFrame(xk)) + vk(12 * i - 5 : 12 * i);
end


function xk1 = f(xk, wk)
xk1 = xk + wk;

