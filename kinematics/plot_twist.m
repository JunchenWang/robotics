function plot_twist(v, name)
% plot twist
if nargin < 2
    name = [];
end
sm = twist2sm(v);
hold on;
axis equal;
for i = 1 : size(v, 1)
    if isempty(name)
        dude = strcat('\leftarrow ','twist', num2str(i));
    elseif iscell(name)
        dude = strcat('\leftarrow ',name{i});
    else
        dude =  strcat('\leftarrow ',name);
    end
    scale = sm(i, 8);
    if isinf(sm(i, 7)) % pure translation
        sty = 'r-';
    elseif  sm(i, 7) == 0 % pure rotation
        sty = 'g-';
    else
        sty = 'black-'; % screw motion
    end
    if scale <= eps
        scale = 1;
        sty = strcat(sty,'-');
    end
    quiver3(sm(i, 1), sm(i, 2), sm(i, 3), ...
            sm(i, 4), sm(i, 5), sm(i, 6), scale, sty, 'LineWidth', 2);
     text(sm(i, 1) + scale * sm(i, 4) / 2, sm(i, 2) + scale * sm(i, 5) / 2, ...
         sm(i, 3) + scale * sm(i, 6) / 2, dude);
end
hold off;