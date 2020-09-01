function f = ConcatenateFrame(f2, f1)
r2 = f2(1 : 3);
t2 = f2(4 : 6);
r1 = f1(1 : 3);
t1 = f1(4 : 6);

r = Concatenate_r(r2, r1);
t = ApplyRotationToPoint(r2, t1) + t2;
f = [r;t];
