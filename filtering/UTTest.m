function UTTest
f = @(x) [x(1) * cos(x(2)); x(1) * sin(x(2))];
J = @(x) [cos(x(2)) -x(1) * sin(x(2)); sin(x(2)) x(1) * cos(x(2))];
r0 = 1;
theta0 = pi / 2;
x0 = [r0; theta0];
Sigma = [0.02^2, 0; 0, (pi / 180 * 15)^2];

N = 100000;
PolarPts = repmat(x0, 1, N) + chol(Sigma) * randn(2, N);
EulerPts = zeros(2, N);
for i = 1 : N
    EulerPts(:, i) = f(PolarPts(:, i));
end
E_sigma_mon = cov(EulerPts');
E_mu_mon = mean(EulerPts, 2);

L_sigma = J(x0) * Sigma * J(x0)';
L_mu = f(x0);
plot(EulerPts(1, :), EulerPts(2, :), 'r+');
hold on;
DrawUncertainty(L_mu, L_sigma, 'b');
DrawUncertainty(E_mu_mon, E_sigma_mon, 'y');

% UT
[ut_mu ut_sigma] = UT(x0, Sigma, f);
DrawUncertainty(ut_mu, ut_sigma, 'g');