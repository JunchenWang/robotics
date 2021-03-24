function T = forward_kin_dh(dh_table)
% dh_table: n x 4 matrix
n = size(dh_table,1);
T = eye(4);
for i = 1 : n
    T = T * rttr(dh_table(i,:));
end