lbr = importrobot('urdf\iiwa7\iiwa7.urdf');
lowers = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
uppers = -lowers;
framerate = 25;
r = rateControl(framerate);
% angles = [20, 0, 40, 50, 60, 70, 80] / 180 * pi;
% T = forward_kin_kuka(angles);
% pose = [
% -0.016776 0.275637 0.961115
% 0.004811 0.961262 -0.275595
% -0.999848 -0.000000 -0.017452
% 779.442119 -223.501431 415.832349
% ];
R = T(1:3,1:3);
t = T(1:3,4);
[angles, bounds] = inverse_kin_kuka(R, t, cfg,lowers, uppers);
q0 = homeConfiguration(lbr);
for i = 1 : 7
    q0(i).JointPosition=angles(i);
end
show(lbr,q0);

% angles = inverse_kin_kuka_kesai(R, t, cfg, 0);
n = length(bounds) / 2;
m = 1000;
cnt = 1;
ang = zeros(m,7,n);

    figure;
    hold on;
    grid on;
    xlabel 'arm angle \psi';
    ylabel 'joint angle \theta'
    psi = linspace(-pi, pi, m);
for kesai = psi
angles = inverse_kin_kuka_kesai(R, t, cfg, kesai);
    ang(cnt,:, k)=angles;
    cnt = cnt + 1;
end
for i = 1 : 7
plot(psi/pi*180, ang(:,i,k)/pi*180,'linewidth',2, 'color', color(i));
end
legend('\theta_1', '\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7');
for k = 1 : n
    plot([bounds(2*k-1),bounds(2*k-1)]/pi*180,[-180,180],'b-.','linewidth',1);
    plot([bounds(2*k),bounds(2*k)]/pi*180,[-180,180],'b-.','linewidth',1);
end
hold off;

hold on;
color='cmkrgby';
for k = 1 : n
    cnt = 1;
    psi = linspace(bounds(2*k-1),bounds(2*k),m);
    figure;
    hold on;
    grid on;
    xlabel 'arm angle \psi';
    ylabel 'joint angle \theta'
for kesai = psi
    angles = inverse_kin_kuka_kesai(R, t, cfg, kesai);
    ang(cnt,:, k)=angles;
    cnt = cnt + 1;
%     for i = 1 : 7
%         q0(i).JointPosition=angles(i);
%     end
%     show(lbr,q0);
%     hold on;
%     waitfor(r);
end
for i = 1 : 7
plot(psi/pi*180, ang(:,i,k)/pi*180,'linewidth',2, 'color', color(i));
end
legend('\theta_1', '\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7');
hold off;
end