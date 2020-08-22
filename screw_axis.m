function sa = screw_axis(p, dir, h)
% p and dir are nx3, representing a point and a direction
% direction does not need to be normalized
if iscolumn(p)
    p = p';
    dir = dir';
end
n = size(p,1);
sm = cat(2, p, dir, ones(n, 1) * h, ones(n, 1));
sa = sm2twist(sm);