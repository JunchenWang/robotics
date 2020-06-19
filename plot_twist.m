function plot_twist(v)
sm = twist2sm(v);
hold on;
axis equal;
quiver3(sm(:, 1), sm(:, 2), sm(:, 3), sm(:, 4), sm(:, 5), sm(:, 6), sm(:, 7), 'black-', 'LineWidth', 2);
hold off;