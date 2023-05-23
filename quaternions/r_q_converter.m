function qf = r_q_converter(fq)
% one to one mapping in the range [0 ~ pi]
fq = fq(:);
if length(fq) == 3
    theta = norm(fq);
    if theta > 0
        qf = [cos(theta / 2); sin(theta / 2) * fq / theta];
    else
        qf = [1, 0, 0, 0]';
    end
    return;
else
    if fq(1) < 0
        fq = -fq;
    end
    if (fq(1) >= 1)
        qf = [0;0;0];
    else
        theta = acos(fq(1));
        qf = 2 * theta * fq(2:4) / sin(theta);
    end
end
