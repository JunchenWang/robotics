function plot_frame(T, scale)
if nargin < 1
    T = make_tform(eye(3), [0 0 0]);
end
if nargin < 2
    scale = 1;
end
n = size(T, 3);
x = reshape(T(1:3,1,:), 3, n)';
y = reshape(T(1:3,2,:), 3, n)';
z = reshape(T(1:3,3,:), 3, n)';
o = reshape(T(1:3,4,:), 3, n)';
hold on;
axis equal;
quiver3(o(:, 1), o(:, 2), o(:, 3), x(:, 1), x(:, 2), x(:, 3), scale, 'r-', 'LineWidth', 2);
quiver3(o(:, 1), o(:, 2), o(:, 3), y(:, 1), y(:, 2), y(:, 3), scale, 'g-', 'LineWidth', 2);
quiver3(o(:, 1), o(:, 2), o(:, 3), z(:, 1), z(:, 2), z(:, 3), scale, 'b-', 'LineWidth', 2);
hold off;