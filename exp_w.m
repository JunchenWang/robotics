function R = exp_w(w)
if iscolumn(w)
    w = w';
end
n = size(w, 1);
R = zeros(3, 3, n);
for i = 1 : n
    theta = norm(w(i, :));
    if theta <= eps
        R(:,:,i) = eye(3);
    else
        w_hat = w(i,:) / theta;
        so_w_hat = so_w(w_hat);
        so_w_hat2 = so_w_hat * so_w_hat;
        R(:,:,i) = eye(3) + sin(theta) * so_w_hat + (1 - cos(theta)) * so_w_hat2;
    end
end