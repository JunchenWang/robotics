function X = damping_least_square(A, B, lambda, thr)
% min ||AX - B||^2 + lambda ||X||^2
AAt = A * A';
if cond(AAt) > thr
    X = (A' * A + lambda^2 * eye(size(A, 2))) \ (A' * B);
else
    X = A' * (AAt \ B);
end