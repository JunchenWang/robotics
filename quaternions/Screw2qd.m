function qd = Screw2qd(theta, d, l, m)
l = l(:);
m = m(:);
theta_2 = theta / 2;
c = cos(theta_2);
s = sin(theta_2);
qd = [c, -d / 2 * s
      s * l, s * m + d / 2 * c * l];
      