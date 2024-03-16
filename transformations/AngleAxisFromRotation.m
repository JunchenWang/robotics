function [r, theta] = AngleAxisFromRotation(R)
I = eye(3);
A = R - I;
[~, ~, V] = svd(A);
v = V(:, end);
v_hat = [R(3, 2) - R(2, 3), R(1, 3) - R(3, 1), R(2, 1) - R(1, 2)]';
phi = atan2(dot(v, v_hat), trace(R) - 1);
r = phi * v;
if nargout == 2
    theta = norm(r);
    if theta > 0
        r = r / theta;
    else
        r = [0, 0, 1]';
    end
end