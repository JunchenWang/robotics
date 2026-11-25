m = 1;
M = 1;
l = 0.5;
g = 9.81;
A = [0, 0, 1, 0; 
     0, 0, 0, 1; 
     0, 3 * m * g / (4 * M + m), 0, 0;
     0, 6 * (M + m) * g / (l * (4 * M + m)), 0, 0];
B = [0;0;4/ (4 * M + m); 6 / (l * (4 * M + m))];
C = [1, 0, 0, 0];
D = 0;

p = cal_second_order_pole(0.0001, 1);
p = [4 * real(p(1)), 3.5 * real(p(1)), p];
K = place(A, B, p);
%% 极点配置
% L = place(A', C', p)';
kf = -1 / (C*((A - B*K) \ B));
sys_d = c2d(ss(A,B,C,D), 0.002, 'zoh');