function J = DiffJacobian(f, x)
n = length(x);
F = f(x);
m = length(F);
J = zeros(m, n);

for i = 1 : n
    xx = x;
    epsilon = 1e-9 * abs(x(i));
    xx(i) = x(i) + epsilon;
    FF = f(xx);
    J(:,i) = (FF - F) / epsilon;
end