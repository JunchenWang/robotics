function TestKalmanFilter
%Test Linear Kalman Filter
A = eye(4);
A(1, 3) = 1;
A(2, 4) = 1;
H = [1 0 0 0; 0 1 0 0];
Q = diag([0.02, 0.02, 0.01, 0.005].^2);
R = diag([0.2, 0.2].^2);
x0 = [10, 10, 1, 0]';
T = 100;
x = zeros(4, T + 1);
x(:, 1) = x0;
y = zeros(2, T);
w = sqrt(Q) * randn(4, T);
v = sqrt(R) * randn(2, T);
for i = 2 : T + 1
    x(:, i) = A * x(:, i - 1) + w(:, i - 1); 
    y(:, i - 1) = H * x(:, i) + v(:, i - 1);
end
x0 = [10 10 0.8 0.1]';
P0 = eye(4);
P0(3, 3) = 0.2;
P0(4, 4) = 0.1;
figure;
plot(x(1, :), x(2, :), 'rs-');
hold on;
plot(y(1, :), y(2, :), 'b*-');
[x P] = KalmanFilter(A, H, x0, P0, y, R, Q);
plot(x(1, :), x(2, :), 'g<-');
% DrawUncertainty(x, P, [], 2, [], 2);



% DC simulator
Q = 1e-5;
R = 0.1^2;
interval = 0.05;
N = 200;
Y = 0.5;
y = Y + sqrt(R) * randn(1,N);
t = interval * [1 : N];
A = 1;
H = 1;
x0 = 0.3;
P0 = 0.2;
[x P] = KalmanFilter(A, H, x0, P0, y, R, Q);
figure;
plot(t, y, 'r');
hold on;
plot([0 t], Y * ones(1, N + 1), 'm');
% DrawUncertainty(x, P, 'c', 1, interval);
plot([0 t], x, 'g*');

    
    
    
