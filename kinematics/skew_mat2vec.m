function x = skew_mat2vec(X)
% x: m by 3 matrix
% X: 3 by 3 by m matrix
m = size(X, 3);
x = zeros(m, 3);
for i = 1 : m
    x(i, :) = [-X(2, 3, i) X(1, 3, i) -X(1, 2, i)];
end

