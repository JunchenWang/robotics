function T =UR_forward_kin(theta) %theta 单位为rad
% deprecated, do not use
dh_table=[0,0,0.1625,theta(1);
          pi/2,0,0,pi+theta(2);
          0,0.425,0,theta(3);
          0,0.3922,0.1333,pi+theta(4);
          pi/2,0,0.0997,theta(5);
          -pi/2,0,0.0996,theta(6)]; 
 T=forward_kin_dh(dh_table) ;
