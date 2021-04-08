% function trajectory_test
aviobj=VideoWriter('example');%新建叫example.avi的文件
open(aviobj); %打开example.avi的文件
framerate = 30;
r = rateControl(framerate);
lbr = importrobot('..\urdf\iiwa7\iiwa7.urdf');
lowers = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
uppers = -lowers;
R0 = [0, 0, 1;0 1 0;-1 0 0];
t0 = [705.75, -300, 365.31]';
cfg = [1, -1, -1];
numSamples = 100;
pts = zeros(3,numSamples);
q0 = homeConfiguration(lbr);
color='cmkrgby';
show(lbr,q0);
axis([-0.3, 1, -0.5, 0.5, 0, 1.35]);
campos([7.84, 8.17,2.05]); camva(6.9); camtarget([0.16, -0.18, 0.55]);
[angles, bounds,kesai] = inverse_kin_kuka(R0, t0, cfg,lowers, uppers);
[q,qd,qdd,tSamples,pp] = trapveltraj([0, 0, 0, 0, 0, 0, 0; angles]',numSamples);
for k = 1 : numSamples
    
    for i = 1 : 7
        q0(i).JointPosition=q(i,k);
    end
    T = forward_kin_kuka(q(:,k));
    pts(:,k) = T(1:3,4);
    show(lbr,q0);
    axis([-0.3, 1, -0.5, 0.5, 0, 1.35]);
    campos([7.84, 8.17,2.05]); camva(6.9); camtarget([0.16, -0.18, 0.55]);
    hold on;
    plot3(pts(1,1:k)/1000,pts(2,1:k)/1000,pts(3,1:k)/1000, 'LineWidth',2, 'Color','r');
    hold off;
%     waitfor(r);
    currFrame = getframe;
    writeVideo(aviobj,currFrame);
end
pts2 = pts;
% numSamples = 500;
pts = zeros(3,numSamples);
y = linspace(0, 600, numSamples);
z = 100*cos(2*pi/300.*y)-100;
ang = zeros(numSamples, 7);
% kesai = 0;
for k = 1 : numSamples
    t = t0 + [0;y(k);z(k)];
    [angles, bounds, kesai_] = inverse_kin_kuka(R0, t, cfg,lowers, uppers);
    if isInRange(bounds, kesai)
        angles = inverse_kin_kuka_kesai(R0, t, cfg,kesai);
    else
        kesai = kesai_;
    end
    ang(k,:) = angles;
    T = forward_kin_kuka(angles);
    pts(:,k) = T(1:3,4);
    for i = 1 : 7
        q0(i).JointPosition=angles(i);
    end
    show(lbr,q0);
    axis([-0.3, 1, -0.5, 0.5, 0, 1.35]);
    campos([7.84, 8.17,2.05]); camva(6.9); camtarget([0.16, -0.18, 0.55]);
    hold on;
    plot3(pts2(1,:)/1000,pts2(2,:)/1000,pts2(3,:)/1000, 'LineWidth',2, 'Color','r');
    plot3(pts(1,1:k)/1000,pts(2,1:k)/1000,pts(3,1:k)/1000, 'LineWidth',2,'Color','b');
    hold off;
%     waitfor(r);
    currFrame = getframe;
    writeVideo(aviobj,currFrame);
end
close(aviobj);
hold on;
q0 = homeConfiguration(lbr);
show(lbr,q0);
hold off;
figure;
hold on;
xlabel 'sample point';
ylabel 'joint angle \theta'
for i = 1 : 7
    plot(1:numSamples, ang(:,i)/pi*180,'linewidth',2, 'color', color(i));
end
legend('\theta_1', '\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7');
grid on;
hold off;
figure;
xyz = [705.75*ones(1,numSamples);t0(2)+y;t0(3)+z];
plot3(xyz(1,:),xyz(2,:),xyz(3,:),'r',pts(1,:),pts(2,:),pts(3,:),'b');
err = sqrt(sum(sum((xyz-pts).^2,1))/numSamples);
