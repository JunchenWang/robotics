function finv = InverseFrame(f)
r = f(1 : 3);
t = f(4 : 6);
R = RotationByAxisAngleRep(r);
t = -R' * t;
finv = zeros(6, 1);
finv(1 : 3) = -r;
finv(4 : 6) = t;