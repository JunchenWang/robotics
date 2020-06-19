function invT = tform_inv(T)
% T: 4x4xn
n = size(T, 3);
invT = zeros(4, 4, n);
for i = 1 : n
    invT(:, :, i) = [T(1:3, 1:3, i)', -T(1:3, 1:3, i)' * T(1:3, 4, i); 0 0 0 1];
end