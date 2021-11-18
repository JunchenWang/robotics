N = size(data,1);
rcm = x(4:6)';
dir = x(7:9)';
dir = dir / norm(dir);
for i = 1 :50: N
    R = EulZYX2R(data(i, 1:3));
    t = data(i, 4:6)';
    p = R*rcm + t;
    d = R *dir;
    p = p - 600 * d;
    quiver3(p(1), p(2), p(3), d(1),d(2),d(3),500);
    hold on;
end
axis equal