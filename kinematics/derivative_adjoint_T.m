function [dAdT, AdT] = derivative_adjoint_T(T, dT)
R = T(1:3,1:3);
t = T(1:3, 4);
dR = dT(1:3,1:3);
dt = dT(1:3, 4);
dAdT = [dR, zeros(3); so_w(dt)*R+so_w(t)*dR, dR];
AdT = [R, zeros(3); so_w(t) * R, R];    
    