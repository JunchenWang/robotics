lbr = importrobot('urdf\iiwa7\iiwa7.urdf');
% lowers = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
% uppers = -lowers;
framerate = 25;
r = rateControl(framerate);
R = T(1:3,1:3);
t = T(1:3,4);
[angles, bounds] = inverse_kin_kuka(R, t, cfg,lowers, uppers);
q0 = homeConfiguration(lbr);
for i = 1 : 7
    q0(i).JointPosition=angles(i);
end
show(lbr,q0);

n = length(bounds) / 2;
m = 10;
ang = zeros(m,7,n);
color='cmkrgby';
for k = 1 : n
    cnt = 1;
    psi = linspace(bounds(2*k-1),bounds(2*k),m);
    for kesai = psi
        angles = inverse_kin_kuka_kesai(R, t, cfg, kesai);
        ang(cnt,:, k)=angles;
        cnt = cnt + 1;
        for i = 1 : 7
            q0(i).JointPosition=angles(i);
        end
        show(lbr,q0);
        hold on;
        waitfor(r);
    end
end