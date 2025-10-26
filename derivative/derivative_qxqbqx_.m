function [Jx Jb] = derivative_qxqbqx_(qx, qb)
% qx * qb * qx_
ax = qx(1);
vx = qx(2:4);
ab = qb(1);
vb = qb(2:4);
vbvx = dot(vb, vx);
J1 = [vx(1) * vb'; vx(2) * vb'; vx(3) * vb'] + vbvx * eye(3);
J2 = -2 * SkewMatrix(ax * vb);
J3 = -derivative_AxCrossBx_x(SkewMatrix(vb)', eye(3), vx);
Jx = [2 * ax * ab, 2 * ab * vx'
     2 * ax * vb - 2 * cross(vb, vx), J1 + J2 + J3];
 qxqb = ConcatenateQuaternions(qx, qb);
 J1 = derivative_q2starq1(qxqb, Conjugate_q(qx));
 [J2 J3] = derivative_q2starq1(qx, qb);
 Jb = J1 * J3;
    
   