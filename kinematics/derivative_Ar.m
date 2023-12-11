function dA = derivative_Ar(r, dr)
% r = 0处的导数如何定义：
% dr = rand(3,1); r = 1e-8 * (dr + [1,0,0]'); derivative_Ar(r, dr)+ SkewMatrix(dr)
nr = norm(r);
if nr == 0
    dA = -so_w(dr);
    return;
end
dnr = dot(dr,r) / nr;
skr = so_w(r);
sk2r = skr * skr;
dskr = so_w(dr);
dsk2r = so_w(dr) * so_w(r) + so_w(r) * so_w(dr);
a1 = (cos(nr) - 1) / nr^2;
a2 = (nr - sin(nr)) / nr^3;
da1 = -sin(nr) * dnr / nr^2 + 2 * (1 - cos(nr)) / nr^3 * dnr;
da2 = (dnr - cos(nr) * dnr) / nr^3 - 3 * (nr - sin(nr)) / nr^4 * dnr;
dA = a1 * dskr + da1 * skr + a2 * dsk2r + da2 * sk2r;