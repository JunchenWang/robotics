function T = forward_kin_kuka(angles)
% angles: A1-An
d1 = 340;
d3 = 400;
d5 = 400;
d7 = 126;
dh_table = [0 0 d1 0; 
            -pi/2 0 0 0;
            pi/2 0 d3 0;
            pi /2 0 0 0;
            -pi/2 0 d5 0;
            -pi/2 0 0 0;
            pi/2 0 d7 0];
n = length(angles);
T = eye(4);
for i = 1 : n
    dh_table(i,4) = angles(i);
    T = T * rttr(dh_table(i,:));
end