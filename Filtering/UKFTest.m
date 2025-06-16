function UKFTest
%reference: unscented Filtering and nonlinear Estimation
x0 = [6500.4; 349.14;-1.8093;-6.7967;0.6932];
MC = 1;
T = 200;
obs_freq = 10;
N = T * obs_freq;
x = zeros(5 * MC, N + 1);
x(:, 1) = repmat(x0, MC, 1);
y = zeros(2 * MC, N);
% P0 = 1e-6 * eye(5);
% P0(5,5) = 0;
Q = 2.4064e-5 * eye(3);
Q(3, 3) = 0;
R = [1, 0; 0, 27e-3];

for j = 1 : MC
    v = sqrt(R) * randn(2, N);
    w = sqrt(Q) * randn(3, N);
    for i = 2 : N + 1
        x(5 * (j - 1) + 1 : 5 * j, i) = kalman_f(x(5 * (j - 1) + 1 : 5 * j, i - 1), w(:, i - 1));
        y(2 * (j - 1) + 1 : 2 * j, i - 1) = kalman_h(x(5 * (j - 1) + 1 : 5 * j, i), v(:, i - 1));
    end
end

plot(mean(x(1 : 5 : 5 * MC, 1 : obs_freq : end), 1), mean(x(2 : 5 : 5 * MC, 1 : obs_freq : end), 1), 'b<');
hold on;
 plot(6000 + cos(y(2, 1 : obs_freq : end)) .* y(1, 1 : obs_freq : end),  sin(y(2, 1 : obs_freq : end)) .* y(1, 1 : obs_freq : end), 'y*');
%  axis([6300 6510 -200 500]);
% x0 = [6500.4; 349.14;-1.8093;-6.7967;0];
hold on;
P0 = 1e-6 * eye(5);
% P0(5, 5) = 2;
% Q = 2.4064e-5 * eye(3);
Q(3, 3) = 1e-6;
observeIndex = 1;
[x P] = UKalmanFilter(@kalman_f, @kalman_h, x0, P0, y(2 * (observeIndex - 1) + 1 : 2 * observeIndex, :), R, Q);
plot(x(1, 1 : obs_freq : end), x(2, 1 : obs_freq : end), 'r*');
[x P] = EKalmanFilter(@kalman_f, @kalman_h, x0, P0, y(2 * (observeIndex - 1) + 1 : 2 * observeIndex, :), R, Q);
plot(x(1, 1 : obs_freq : end), x(2, 1 : obs_freq : end), 'mo');
% plot(t, y, 'r');
% hold on;
% plot([x0 t], Y * ones(1, N + 1), 'm');
% plot([0 t], [x0 x], 'b');
% DrawUncertainty(x, P, 'c', 1.5, interval);
% DrawUncertainty([0;x0], diag([P0, P0]), 'y');
% hold off;

function [xk1 JA JW] = kalman_f(xk, w)
%reference: unscented Filtering and nonlinear Estimation
xk = xk(:);
R = sqrt(xk(1).^2 + xk(2).^2);
V = sqrt(xk(3).^2 + xk(4).^2);
beta0 = -0.59783;
H0 = 13.406;
Gm0 = 3.9860e5;
R0 = 6374;
beta = beta0 * exp(xk(5));
G = -Gm0 / R.^3;
D = beta * exp((R0 - R) / H0) * V;
deta = [xk(3); xk(4); D * xk(3) + G * xk(1) + w(1); D * xk(4) + G * xk(2) + w(2); w(3)];
xk1 = xk + deta;
JR = [xk(1) xk(2) 0 0 0] / R;
JV = [0 0 xk(3) xk(4) 0] / V;
Jbeta = beta * [0 0 0 0 1];
JG = 3 * Gm0 / R.^4 * JR;
JD = Jbeta * exp((R0 - R) / H0) * V... 
     + beta * exp((R0 - R) / H0) * (-JR) / H0 * V ...
     + beta * exp((R0 - R) / H0) * JV;
Jdeta1 = [0 0 1 0 0];
Jdeta2 = [0 0 0 1 0];
Jdeta3 = JD * xk(3) + D * [0 0 1 0 0] + JG * xk(1) + G * [1 0 0 0 0];
Jdeta4 = JD * xk(4) + D * [0 0 0 1 0] + JG * xk(2) + G * [0 1 0 0 0];
Jdeta5 = [0 0 0 0 0];
JA = eye(5) + [Jdeta1; Jdeta2; Jdeta3; Jdeta4; Jdeta5];
JW = [0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1];

function [yk JX JV] = kalman_h(xk, v)
%reference: unscented Filtering and nonlinear Estimation
xr = 6000;
yr = 0;
xk = xk(:);
xx = xk(1) - xr;
yy = xk(2) - yr;
r = sqrt(xx.^2 + yy.^2);
theta = atan(yy / xx);
yk = [r + v(1); theta + v(2)];
JX1 = 1 / r * [xx, yy, 0, 0, 0];
JX2 = 1 / (1 + (yy / xx).^2) * [-yy / xx.^2, 1 / xx, 0, 0, 0];
JX = [JX1; JX2];
JV = eye(length(v));
