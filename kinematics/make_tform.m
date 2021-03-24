function T = make_tform(R, t)
% R: 3x3xn
% t: nx3
if nargin == 0
    T = [eye(3), zeros(3, 1); 0 0 0 1];
    return;
end
n = size(R, 3);
T = zeros(4, 4, n);
for i = 1 : n
    T(:, :, i) = [R(:, :, i), t(i, :)'; 0 0 0 1];
end