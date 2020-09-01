function [R t FRE_2 s] = OrientationProblem(x, y)
% y = R*x + t
% y,x: 3 by N matrix, N pionts
nFrames = size(x, 2);
x0 = sum(x, 2) / nFrames;
y0 = sum(y, 2) / nFrames;
x_ = x - repmat(x0, 1, nFrames);
y_ = y - repmat(y0, 1, nFrames);
M = y_ * x_';
s = norm(y_, 'fro') / norm(x_, 'fro');
[U D V] = svd(M);
R = U * diag([1,1,det(U*V')]) * V';
if nargout == 4
    t = y0 - s * R * x0;
    FRE_2 = sum(sum((s * R * x + repmat(t, 1, nFrames) - y).^2)) / nFrames;
else 
    t = y0 - R * x0;
    FRE_2 = sum(sum((R * x + repmat(t, 1, nFrames) - y).^2)) / nFrames;
end