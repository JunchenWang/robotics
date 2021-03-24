function T = rttr(dh_row)
% dh_row: alpha, a, d, theta
alpha = dh_row(1);
a = dh_row(2);
d = dh_row(3);
theta = dh_row(4);
ct = cos(theta);
st = sin(theta);
cpha = cos(alpha);
spha = sin(alpha);
T = [ct, -st, 0, a; 
    st * cpha, ct * cpha, -spha, -d * spha; 
    st * spha, ct * spha, cpha, d * cpha;
    0, 0, 0, 1];

