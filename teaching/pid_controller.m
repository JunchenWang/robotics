function u = pid_controller(t, y, yd, Kp, Ki, Kd)
persistent sum;
persistent pret;
if isempty(sum)
    sum = 0;
    pret = 0;
end
ydt = yd(t);
e = ydt(1) - y(1);
ed = ydt(2) - y(2);
sum = sum + e * (t - pret);
u = Kp * e + Kd * ed + Ki * sum;
pret = t;

