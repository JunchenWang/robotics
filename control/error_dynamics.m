m = 1;
b = 10;
k = 25;
y0 = [1, 0]';
T = [0, 10];
w = sqrt(k / m);
xi = b / (2*sqrt(k*m));
[t,y] = ode45(@(t, y) err(t, y, m, k, b),T, y0);
% plot(t, y(:,1), t, y(:,2));
figure;
plot(t, y(:,1));
function dydt = err(t,y,m, k, b)
if m ~= 0
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = (-b * y(2) - k * y(1)) / m;
else
    dydt = -k / b * y;
end
end