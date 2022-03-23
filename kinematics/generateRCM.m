function T = generateRCM(d,a,b,c,p_rcm)
x_rcm = [0, 1,0]';
z_rcm = [0, 0, 1]';
y_rcm = cross(z_rcm, x_rcm);
 
a = deg2rad(a);
b = deg2rad(b);
c = deg2rad(c);

St = sm2twist([p_rcm', z_rcm', inf, 1]) * d;
% p_rcm = p_rcm - d * z_rcm;
Sz = sm2twist([p_rcm', z_rcm', 0, 1]) * a;
Sx = sm2twist([p_rcm', x_rcm', 0, 1]) * b;
Sy = sm2twist([p_rcm', y_rcm', 0, 1]) * c;
T = exp_twist(Sx)*exp_twist(Sy)*exp_twist(Sz)*exp_twist(St);
