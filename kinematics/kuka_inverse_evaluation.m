lbr = importrobot('..\urdf\iiwa7\iiwa7.urdf');
lowers = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
uppers = -lowers;
% R = [
%    -0.3585    0.4969    0.7903
%    -0.8523   -0.5196   -0.0599
%     0.3809   -0.6951    0.6098];
% t = [779.442119, -223.501431, 415.832349]';
cfg = [1,1,1];
% R = T(1:3,1:3);
% t = T(1:3,4);
[angs, bounds] = inverse_kin_kuka(R, t, cfg,lowers, uppers);
q0 = homeConfiguration(lbr);
color='cmkrgby';
for i = 1 : 7
    q0(i).JointPosition=angs(i);
end
show(lbr,q0);

n = length(bounds) / 2;
m = 1000;
cnt = 1;
ang = zeros(m,7);

figure;
hold on;
grid on;
xlabel 'arm angle \psi';
ylabel 'joint angle \theta'
psi = linspace(-pi, pi, m);
for kesai = psi
    angs = inverse_kin_kuka_kesai(R, t, cfg, kesai);
    ang(cnt,:)=angs;
    cnt = cnt + 1;
end
for i = 1 : 7
    plot(psi/pi*180, ang(:,i)/pi*180,'linewidth',2, 'color', color(i));
end
legend('\theta_1', '\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7');
for k = 1 : n
    plot([bounds(2*k-1),bounds(2*k-1)]/pi*180,[-180,180],'b-.','linewidth',1);
    plot([bounds(2*k),bounds(2*k)]/pi*180,[-180,180],'b-.','linewidth',1);
end
hold off;
