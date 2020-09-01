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
     if (fq(1) == 1)
         qf = [0;0;0];
     else
        theta = acos(fq(1));
         qf = 2 * theta * fq(2:4) / sin(theta);
     end

%      v = fq(2 : 4);
%      if norm(v) == 0 
%          qf = [0 0 0]';
%      else       
%          qf = 2 * sign(fq(1)) * asin(norm(v)) / norm(v) * v;
%      end
%     theta = 2 * atan2(norm(fq(2 : 4)), fq(1));
%     if theta == 0
%         qf = [0 0 0]';
%     else
%         qf = fq(2 : 4) / sin(theta / 2) * theta;
%     end
%       v = fq(2 : 4);
%       qf = 2 * sign(fq(1)) * asin(norm(v)) / norm(v) * v;
%     if theta > pi
%         qf = -qf;
%         theta = 2 * pi - theta;
%     end
%     qf = qf * theta;
end

function f = sign(x)
if x >= 0
    f = 1;
else
    f = -1;
end