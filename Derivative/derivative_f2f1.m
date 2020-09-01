function [Jf2 Jf1] = derivative_f2f1(f2, f1)
f2 = f2(:);
f1 = f1(:);
r2 = f2(1 : 3);
r1 = f1(1 : 3);
t1 = f1(4 : 6);

[Jr2 Jr1] = derivative_r2r1(r2, r1);

Jf2 = [Jr2 zeros(3); derivative_rx_r(r2, t1) eye(3)];
Jf1 = [Jr1 zeros(3); zeros(3) RotationByAxisAngleRep(r2)];