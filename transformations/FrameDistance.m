function d = FrameDistance(f1, f2)
n = size(f1, 2);
d = zeros(2, n);
for i = 1 : n
    [dr dt] = MatrixDistanceSquare(f1(:, i), f2(:, i));
    d(:, i) = [dr; dt];
end
d = sqrt(mean(d, 2));
