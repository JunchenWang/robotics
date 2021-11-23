function u = pid_controller(t, y, yd, Kp, Ki, Kd)
persistent sum;
persistent pret;
if isempty(sum)
    sum = 0;
    pret = 0;
end
e = yd(1) - y(1);
ed = yd(2) - y(2);
sum = sum + e * (t - pret);
u = Kp * e + Kd * ed + Ki * sum;
pret = t;

