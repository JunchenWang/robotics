function Adv = adjoint_twist(v)
if iscolumn(v)
    v = v';
end
n = size(v, 1);
Adv = zeros(6, 6, n);
for i = 1 : n
    Adv(1:3, 1:3, i) = so_w(v(i, 1:3));
    Adv(4:6, 4:6, i) = Adv(1:3, 1:3, i);
    Adv(4:6, 1:3, i) = so_w(v(i, 4:6));
end