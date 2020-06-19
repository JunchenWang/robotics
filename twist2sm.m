function sm = twist2sm(v)
% sm: q (point), s (u-axis) , theta, h (pitch)
% inf h indicats translation
if iscolumn(v)
    v = v';
end
n = size(v, 1);                                                                                                              
sm = zeros(n, 8);
for i = 1 : n
    theta = norm(v(i, 1:3));           
    if theta <= eps % accuracy issue.
        sm(i, 8) = inf;
        theta = norm(v(i, 4:6));
        if theta <= eps
            warning('v is all zeros');
            sm(i, 6) = 1;
        else
            sm(i,4:7) = [v(i, 4:6) / theta  theta];
        end
    else
        s = v(i, 1:3) / theta;
        h = dot(s, v(i, 4:6)) / theta;
        % q is chosen to be with minimu norm
        q = (pinv(so_w(s)) * (h * s - v(i, 4:6) / theta)')';
        sm(i,:) = [q, s, theta, h];
    end
end
