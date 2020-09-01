function J = derivative_xd_xu(xu, k)
n = size(xu, 2);
J = zeros(2 * n);
for i = 1 : n
    r2 = xu(1, i)^2 + xu(2, i)^2;
    tmp1 = 1 + k(1) * r2 + k(2) * r2^2;
    tmp2 = k(1) + 2 * k(2) * r2;
    J(2 * i - 1 : 2 * i, 2 * i - 1 : 2 * i) = [tmp1 + 2 * xu(1, i)^2 * tmp2 + 2 * k(3) * xu(2, i) + 6 * k(4) * xu(1, i), ...
                                               2 * xu(1, i) * xu(2, i) * tmp2 + 2 * k(3) * xu(1, i) + 2 * k(4) * xu(2, i);
                                               2 * xu(1, i) * xu(2, i) * tmp2 + 2 * k(4) * xu(2, i) + 2 * k(3) * xu(1, i), ...
                                               tmp1 + 2 * xu(2, i)^2 * tmp2 + 2 * k(4) * xu(1, i) + 6 * k(3) * xu(2, i)];
end