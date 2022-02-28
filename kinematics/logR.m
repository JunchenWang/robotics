function w = logR(R)
% rotation R -> axis angle rep
eps = 1e-12;
n = size(R, 3);
w = zeros(n, 3);
for i = 1 : n
    tr = trace(R(:,:,i));
    if abs(tr + 1) <= eps
        w(i,:) = pi / sqrt(2 * (1 + R(1,1,i))) * (R(:,1,i)' + [1 0 0]);
    else
        if (tr - 1) / 2 <= 1 - eps
            theta = acos((tr - 1) / 2);
            v = [R(3,2,:)-R(2,3,:), R(1,3,:)-R(3,1,:), R(2,1,:)-R(1,2,:)] / (2 * sin(theta));
            w(i,:) = theta  * v;
        end
    end
end