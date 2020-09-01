function rn = Normalize_r(r)
r = reshape(r, 3, []);
N = size(r, 2);
rn = zeros(3, N);
for i = 1 : N
    theta = norm(r(:, i));
    v = r(:, i) / theta;
    if (theta > pi)
        v = -v;
        theta = 2 * pi - theta;
        rn(:, i) = theta * v;
    else
        rn(:, i) = r(:, i);
    end
end