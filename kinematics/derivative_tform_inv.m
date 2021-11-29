function [dInvT, invT] = derivative_tform_inv(T, dT)
dR = dT(1:3,1:3);
dt = dT(1:3,4);
R = T(1:3, 1:3);
t = T(1:3, 4);
dInvT = [dR', -dR'*t-R'*dt; zeros(1,4)];
invT = [R', -R'*t; 0 0 0 1];