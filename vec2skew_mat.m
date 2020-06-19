function X = vec2skew_mat(x)
% x: vector or m by 3 matrix
% X: 3 by 3 by m matrix
if isvector(x)
    X = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1), 0];
else
    m = size(x, 1);
    X = zeros(3, 3, m);
    for i = 1 : m
        X(:, :, i) = [0 -x(i, 3) x(i, 2); x(i, 3) 0 -x(i, 1); -x(i, 2) x(i, 1), 0];
    end
end

