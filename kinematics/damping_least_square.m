function X = damping_least_square(A, B, lambda, thr)
% Minimize ||A*X - B||^2 + lambda_sq*||X||^2.
% A positive threshold blends damping in over [0.5*thr, 2*thr].
AAt = A * A';
lambda_sq = lambda^2;
if thr > 0
    if isinf(thr)
        lambda_sq = 0;
    else
        u = min(1, max(0, (cond(AAt) / thr - 0.5) / 1.5));
        weight = u^2 * (3 - 2*u);
        lambda_sq = lambda_sq * weight;
    end
end
X = A' * ((AAt + lambda_sq * eye(size(AAt, 1), 'like', AAt)) \ B);
end
