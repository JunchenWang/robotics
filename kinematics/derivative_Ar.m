function dA = derivative_Ar(r, dr)
nr = norm(r);
if nr < 1e-6
    r = 1e-6 * ones(3,1); % 1e-6 not change, 为了应对r接近于零的奇异位置
    nr = norm(r);
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