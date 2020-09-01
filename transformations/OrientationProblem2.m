function [R t FRE_2 s] = OrientationProblem2(x, y)
N = size(x, 2);
x0 = sum(x, 2) / N;
y0 = sum(y, 2) / N;
x_ = x - repmat(x0, 1, N);
y_ = y - repmat(y0, 1, N);
% M = y_ * x_';
s = norm(y_, 'fro') / norm(x_, 'fro');
CovM = zeros(3);
for i = 1 : N
    CovM = CovM + x(:,i) * y(:,i)';
end
CovM = 1 / N * CovM - x0 * y0';
delta = CovM - CovM';
delta = [delta(2,3); delta(3,1); delta(1,2)];
M = [trace(CovM), delta'; delta, CovM + CovM' - trace(CovM) * eye(3)];
[E D] = eig(M);
D = diag(D);
[d i] = max(D);
q = E(:, i);
R = RMByQuaternion(q);
% [U D V] = svd(M);
% R = U * V';
if nargout == 4
    t = y0 - s * R * x0;
    FRE_2 = sum(sum((s * R * x + repmat(t, 1, N) - y).^2)) / N;
else 
    t = y0 - R * x0;
    FRE_2 = sum(sum((R * x + repmat(t, 1, N) - y).^2)) / N;
end