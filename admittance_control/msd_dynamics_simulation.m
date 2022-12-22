function msd_dynamics_simulation
m = 1;
k = 500;
M = m * eye(3);
K = k * eye(3);
B = 2*sqrt(m * k) * eye(3);
y0 = zeros(6, 1);
y0(4:6) = [0, 0, 1]';
dynamic = @(t, y) msd_dynamics(t, y, M, B, K, @force);
tspan = [0, 2];
[t, y] = ode45(dynamic, tspan, y0);
plot(t, y(:,1:3)');


function f = force(t)
if t > 0.5 && t < 1.5
    f = [0; 0; 0];
else
    f = [0; 0; 0];
end
function yd = msd_dynamics(t, y, M, B, K, f)
x = y(1: 3);
xd = y(4:6);
xdd = M \ (f(t) - B * xd - K * x);
yd = [xd; xdd];