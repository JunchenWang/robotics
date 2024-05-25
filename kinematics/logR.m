function w = logR(R)
n = size(R, 3);
w = zeros(n, 3);
I = eye(3);
for i = 1 : n
    [~, ~, V] = svd(R(:,:,i) - I);
    v = V(:, end);
    v_hat = [R(3, 2, i) - R(2, 3, i), R(1, 3, i) - R(3, 1, i), R(2, 1, i) - R(1, 2, i)]';
    phi = atan2(dot(v, v_hat), trace(R) - 1);
    w(i,:) = phi * v;
end