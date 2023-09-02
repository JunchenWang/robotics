function G = spatial_inertia_matrix(I, m, com)
% return the spatial matrix G with respect to the frame of I
% I is not necessary defined in center mass
% com is a row vector representing the center mass in I's frame 
G = zeros(6,6);
skcom = so_w(com);
G(1:3,1:3) = I - m * (com * com' * eye(3) - com' * com + skcom * skcom);
G(1:3,4:6) = m * skcom;
G(4:6,1:3) = -G(1:3,4:6);
G(4:6,4:6) = m * eye(3); 