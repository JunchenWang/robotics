function v = logT(T)
n = size(T, 3);
v = zeros(n, 6);
I = eye(3);
for i = 1 : n
    R = T(1:3,1:3,i);
    p = T(1:3,4,i);
    v(i, 1:3) = logR(R);
    theta = norm(v(i, 1:3));
    if theta <= eps
        v(i,:) = [0 0 0 p'];
        continue;
    end
    w = v(i, 1:3) / theta;
    W = so_w(w);
    W2 = W * W;
    v(i, 4:6) = (I - theta / 2 * W + (1 - theta / 2 * cot(theta/2)) * W2) * p;
end