function u = pidVelocityController(t, y, yd, Kp, Ki, Kd)
% yd(1) 期望速度; yd(2) 期望加速度
persistent sum;
persistent pret;
persistent pree;
if isempty(sum)
    sum = 0;
    pret = 0;
    pree = 0;
end
ydt = yd(t);
e = ydt(1) - y(2);
ed = ydt(2) - y(3);
sum = sum + (e + pree) / 2 * (t - pret);%梯形公式计算速度误差积分
u = Kp * e + Kd * ed + Ki * sum;
pret = t;
pree = e;

