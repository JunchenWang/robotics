function logR_test
r = rand(3,1);
r = r / norm(r) * pi; 
R = exp_w(r);
d = 1e-9; % noise
R = R + d;
R = eye(3) + d;
% 真值
disp('真值');
disp(r'); 
% logR
disp('logR');
disp(logR(R)); 
% angleaxis
disp('AngleAxisFromRotation');
disp(AngleAxisFromRotation(R)'); 
% rotm2axangle
disp('rotm2axangle');
ax = rotm2axang(R); 
disp(ax(1:3)*ax(4));
   