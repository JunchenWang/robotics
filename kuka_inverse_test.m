lowers = [-70, -20, -70, -20 ,-70, -20, -75] / 180 * pi;
uppers = -lowers;
N = 100000000;
% error = zeros(N,2);
ang1 = zeros(1000,7);
ang2 = zeros(1000,7);
cnt1 = 0;
cnt2 = 0;
tic;
for i = 1 : N
    angles = round(rand(1,7).*(uppers - lowers - 0.01) + lowers + 0.005,4);
%     if abs(angles(2)) < 1e-1 || abs(angles(6)) < 1e-1
%         continue;
%     end
%     angles(2) = 0;
    T = forward_kin_kuka(angles);
    R = T(1:3,1:3); t = T(1:3,4);
    [angles2, bds] = inverse_kin_kuka(R, t, [sign(angles(2)+eps), sign(angles(4)+eps),sign(angles(6)+eps)], lowers, uppers);
%     angles2 = inverse_kin_kuka_kesai(R, t, [sign(angles(2)+eps), sign(angles(4)+eps),sign(angles(6)+eps)], 0);
    if isempty(bds)
        disp('no kesai');
        cnt1 = cnt1 + 1;
        ang1(cnt1,:) = angles;
        continue;
    end
    T2 = forward_kin_kuka(angles2);
    R2 = T2(1:3,1:3);
    t2 = T2(1:3,4);
    err = [norm(t-t2), norm((R2EulZYX(R) - R2EulZYX(R2))/pi*180)];
    if norm(err) > 0.2
        disp(err);
        cnt2 = cnt2 + 1;
        ang2(cnt2,:) = angles;
    end
end
ang1(cnt1+1:end,:)=[];
ang2(cnt2+1:end,:)=[];
toc;
save;