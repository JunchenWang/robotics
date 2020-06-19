function v = sm2twist(sm)
% sm: q (point), s (u-axis), theta, h (pitch) 
if iscolumn(sm)
    sm = sm';
end
n = size(sm, 1);
v = zeros(n, 6);
for i = 1 : n
    q = sm(i, 1:3);
    s = sm(i, 4:6) / norm(sm(i, 4:6));
    theta = sm(i, 7);
    h = sm(i, 8);
    if isinf(h)
        v(i, :) = [0 0 0 s] * theta;
    else
        v(i, :) = [s, cross(-s, q) + h * s] * theta;
    end
end
