function [s, ds] = trajectory_planner(ds0, t, a)
N = 1000;
tt = linspace(0, t, N);
s = ones(1,N);
ds = ones(1,N);
b = -(a*t+ds0);
c = 1 + ds0^2/(2*a);
delta = b^2 - 4*a*c;
if delta < 0
    error('bad a!');
end
tb = (-b - sqrt(delta)) / (2*a);
v = a*tb;
t1 = (v - ds0) / a;
if t1 < 0
    tb = (1 - ds0^2/(2*a))/(a*t-ds0);
    v = a*tb;
    t1 = (v - ds0) / -a;
    t2 = t - tb;
    for i = 1 : N
        if tt(i) <= t1
            s(i) = -0.5*a*tt(i)^2+ds0*tt(i);
            ds(i) = -a*tt(i)+ds0;
        elseif tt(i) <= t2
            s(i) = -0.5*a*t1^2+ds0*t1 + (tt(i) - t1)*v;
            ds(i) = v;
        elseif tt(i) <= t
            s(i) = 1 - 0.5*a*(t-tt(i))^2;
            ds(i) = a*(t-tt(i));
        end
    end
else
    t2 = t - tb;
    for i = 1 : N
        if tt(i) <= t1
            s(i) = 0.5*a*tt(i)^2+ds0*tt(i);
            ds(i) = a*tt(i)+ds0;
        elseif tt(i) <= t2
            s(i) = 0.5*a*t1^2+ds0*t1 + (tt(i) - t1)*v;
            ds(i) = v;
        elseif tt(i) <= t
            s(i) = 1 - 0.5*a*(t-tt(i))^2;
            ds(i) = a*(t-tt(i));
        end
    end
end
plot(tt, s, tt, ds);