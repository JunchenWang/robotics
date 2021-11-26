function T =UR_forward_kin(theta) %theta 单位为rad
dh_table=[0,0,89.2,theta(1);
          pi/2,0,0,pi+theta(2);
          0,425,0,theta(3);
          0,392,109.3,pi+theta(4);
          pi/2,0,94.75,theta(5);
          -pi/2,0,82.5,theta(6)]; 
 T=forward_kin_dh(dh_table) ;