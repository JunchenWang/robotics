function [theta d l m] = ScrewAxis(R, t)
[l theta] = AngleAxisFromRotation(R);
d = dot(t, l);
c = 0.5 * (t - d * l + 1 / tan(theta / 2) * cross(l, t));
m = cross(c, l);