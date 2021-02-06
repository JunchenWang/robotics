lowers = [-70, -120, -60, -20 ,-30, -80, -175] / 180 * pi;
uppers = -lowers;
N = 100000000;
% error = zeros(N,2);
ang1 = zeros(1000,7);
ang2 = zeros(1000,7);
cnt1 = 0;
cnt2 = 0;
tic;
for i = 1 : N
    angles = rand(1,7).*(uppers - lowers - 0.01) + lowers + 0.005;
    angles(2) = 0;
%     angles=[0.999682065393822,0,0.988769319634779,-5.32450183971401e-06,0.0566354819814592,0.858225343151257,-2.97239768111234];
%     angles=[-0.753650427305537,0,-0.805778530943996,2.43574502084131e-06,-0.175630356359407,-0.351615211427109,-2.23433507516819];
%     angles = [0.713348143969816,0,0.752683006037861,1.00000000000000e-09,0.455933243838771,0.751912728461654,-2.80998666012245];
%     angles(2) = 0;
%     angles(4) = 0.1;
%     angles = [-0.454076870199365,0,-0.127495629757278,1.00000000000000e-09,0.380128637223106,-6.33648735597243e-07,1.25769918465413];
%     if abs(angles(2)) < 1e-1 || abs(angles(6)) < 1e-1
%         continue;
%     end
    T = forward_kin_kuka(angles);
    R = T(1:3,1:3); t = T(1:3,4);
    [angles2, bds] = inverse_kin_kuka(R, t, [sign(angles(2)+eps), sign(angles(4)+eps),sign(angles(6)+eps)], lowers, uppers);
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
toc;
save;