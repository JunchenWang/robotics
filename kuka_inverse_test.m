lower = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
upper = -lower;
N = 100000000;
% error = zeros(N,2);
ang1 = zeros(1000,7);
ang2 = zeros(1000,7);
cnt1 = 0;
cnt2 = 0;
tic;
for i = 1 : N
    angles = round(rand(1,7).*(upper - lower) + lower,3);
%     angles(2) = 0;
%     angles = [0.214611133575105,0.00428898015018175,-2.18653642547971,0.0165379847711349,-1.60252965189329,0.000159335524057891,0.395094142843400];
%     angles = [-1.42500000000000,0,-2.02900000000000,0.795000000000000,-2.44300000000000,-0.866000000000000,0.501000000000000];
    T = forward_kin_kuka(angles);
    R = T(1:3,1:3); t = T(1:3,4);
    [angles2, bds] = inverse_kin_kuka(R, t, [sign(angles(2)+eps), sign(angles(4)+eps),sign(angles(6)+eps)]);
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
    if norm(err) > 1
        disp(err);
        cnt2 = cnt2 + 1;
        ang2(cnt2,:) = angles;
    end
end
ang1(cnt1+1:end,:)=[];
ang2(cnt2+1:end,:)=[];
toc;