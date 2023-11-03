function T = forward_kin_UR(q, param) %theta 单位为rad
% param: |a2|, |a3|, d1, d4, d5, d6
a2 = -param(1);
a3 = -param(2);
d1 = param(3);
d4 = param(4);
d5 = param(5);
d6 = param(6);
dh_table=[0, 0, d1, q(1);
          pi/2, 0, 0, q(2);
          0, a2, 0, q(3);
          0, a3, d4, q(4);
          pi/2, 0, d5, q(5);
          -pi/2, 0, d6, q(6)]; 
 T=forward_kin_dh(dh_table);