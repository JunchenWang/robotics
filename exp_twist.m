function tform = exp_twist(v)
if iscolumn(v)
    v = v';
end
n = size(v, 1);
tform = zeros(4, 4, n);
for i = 1 : n
    tform(4,4,i) = 1;
    sa = screw_axis(v(i, :));
    s = sa(1:6);
    theta = sa(7);
    W = so_w(s(1:3));
    tform(1:3, 1:3, i) = exp_w(v(i, 1:3));
    tform(1:3, 4, i) = (eye(3) * theta + (1 - cos(theta)) * W + (theta - sin(theta)) * W * W) * s(4:6)';
end
