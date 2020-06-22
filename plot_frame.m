function plot_frame(T, name, scale)
% plot frame
if nargin < 1
    T = make_tform(eye(3), [0 0 0]);
end
if nargin < 2
    name = [];
end
if nargin < 3
    scale = 1;
end

n = size(T, 3);
x = reshape(T(1:3,1,:), 3, n)';
y = reshape(T(1:3,2,:), 3, n)';
z = reshape(T(1:3,3,:), 3, n)';
o = reshape(T(1:3,4,:), 3, n)';
hold on;
axis equal;
for i = 1 : n
    if isempty(name)
        dude = strcat('\leftarrow ','Frame', num2str(i));
    elseif iscell(name)
        dude = strcat('\leftarrow ',name{i});
    else
        dude =  strcat('\leftarrow ',name);
    end
    quiver3(o(i, 1), o(i, 2), o(i, 3), x(i, 1), x(i, 2), x(i, 3), scale, 'r-', 'LineWidth', 2);
    quiver3(o(i, 1), o(i, 2), o(i, 3), y(i, 1), y(i, 2), y(i, 3), scale, 'g-', 'LineWidth', 2);
    quiver3(o(i, 1), o(i, 2), o(i, 3), z(i, 1), z(i, 2), z(i, 3), scale, 'b-', 'LineWidth', 2);
    text(o(i, 1), o(i, 2), o(i, 3), dude);
end
hold off;