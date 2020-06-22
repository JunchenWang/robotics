function sa = screw_axis(v)
% v is nx6 each row represents a twist
if iscolumn(v)
    v = v';
end
n = size(v, 1);
sa = zeros(n, 7);
for i = 1 : n
    if norm(v(i,:)) <= eps
        sa(i,6) = 1;
        continue;
    end
    theta = norm(v(i,1:3));
    if theta <= eps
        sa(i, 7) = norm(v(i, 4:6));
        sa(i, 4:6) = v(i, 4:6) / sa(i, 7);
    else
        sa(i, 7) = theta;
        sa(i, 1:6) = v(i,1:6) / theta;
    end
end
