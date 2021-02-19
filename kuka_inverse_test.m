lowers = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
uppers = -lowers;
N = 100000000;
% pts = zeros(3,N);
% error = zeros(N,2);
ang1 = zeros(1000,7);
ang2 = zeros(1000,7);
cnt1 = 0;
cnt2 = 0;
tElapsed = 0;
for i = 1 : N
    angles = rand(1,7).*(uppers - lowers - 0.02) + lowers + 0.01;
    if abs(angles(4)) < 0.1
        continue;
    end
%     angles(4) = 0.1;
    angles(2) = 0;
%     angles = [0, 40, 0, -50, 0, -60, 0] / 180 * pi;
%     angles = [-0.314975724458602,1.73254319472325,1.02805718431992,0.283307704851932,2.70064517360875,-0.256078210760688,1.03034891692248];
%      angles = [   -0.3470, 0,    1.5949,    1.8529 ,  -1.8404, -1.6462    1.6619];%(ap)
% 　　　angles=[-0.886085580026083,1.23567066220902,1.00753288588946,1.85285569167132,-0.364260640510686,-0.502116091080774,1.05779326273372];%ap    
    T = forward_kin_kuka(angles);
    R = T(1:3,1:3); t = T(1:3,4);
%     pts(:,i) = t;
    cfg=[sign(angles(2)+eps), sign(angles(4)+eps),sign(angles(6)+eps)];
    tStart=tic;
    [angles2, bds] = inverse_kin_kuka(R, t, cfg, lowers, uppers);
    tElapsed=tElapsed+toc(tStart);
    if isempty(angles2)
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
save;