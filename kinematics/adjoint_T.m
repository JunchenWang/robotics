function Adt = adjoint_T(tform)
% AdT operator T -> 6x6 mapping
n = size(tform, 3);
Adt = zeros(6, 6, n);
for i = 1 : n
    Adt(1 : 3, 1 : 3, i) = tform(1 : 3, 1 : 3, i);
    Adt(4 : 6, 4 : 6, i) = tform(1 : 3, 1 : 3, i);
    Adt(4 : 6, 1 : 3, i) = so_w(tform(1 : 3, 4, i)) * tform(1 : 3, 1 : 3, i);
end
    
    