function w = logR(R)
eps = 1e-9;
n = size(R, 3);
w = zeros(n, 3);
I = eye(3);
for i = 1 : n
    if norm(R(:,:,i) - I, 'fro') <= eps
        continue;
    end
    tr = trace(R(:,:,i));
    if abs(tr + 1) <= eps
        w(i,:) = pi / sqrt(2 * (1 + R(1,1,i))) * (R(:,1,i)' + [1 0 0]);
    else
        theta = acos((tr - 1) / 2);
        w(i,:) = theta / 2 / sin(theta) * skew_mat2vec(R(:,:,i) - R(:,:,i)');
    end 
end