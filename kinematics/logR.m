function w = logR(R)
% rotation R -> axis angle rep
eps = 1e-7; % do not change
n = size(R, 3);
w = zeros(n, 3);
for i = 1 : n
    tr = trace(R(:,:,i));
    theta = real(acos(complex((tr - 1) / 2)));
    if abs(theta - pi) < eps
        w(i,:) = pi / sqrt(2 * (1 + R(1,1,i))) * (R(:,1,i)' + [1 0 0]);
        % 有时会出现theta = pi，但是上式分母为0的情况
        if any(isnan(w(i,:)))
            [~, ~, V] = svd(R - eye(3));
            w(i,:) = pi * V(:, end);
        end
    elseif theta > eps
        v = [R(3,2,i)-R(2,3,i), R(1,3,i)-R(3,1,i), R(2,1,i)-R(1,2,i)] / (2 * sin(theta));
        normv = norm(v);
        if normv ~= 0
            v = v / normv;
        end
        w(i,:) = theta  * v;
    end
%     if abs(tr + 1) <= eps
%         w(i,:) = pi / sqrt(2 * (1 + R(1,1,i))) * (R(:,1,i)' + [1 0 0]);
%     else
%         if (tr - 1) / 2 <= 1 - eps
%             theta = acos((tr - 1) / 2);
%             v = [R(3,2,:)-R(2,3,:), R(1,3,:)-R(3,1,:), R(2,1,:)-R(1,2,:)] / (2 * sin(theta));
%             w(i,:) = theta  * v;
%         end
%     end
end