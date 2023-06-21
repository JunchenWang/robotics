function dZ = derivative_null_z(J, dJ)
[m, n] = size(J);
dJm = dJ(1:m, 1:m);
dJr = dJ(1:m, m + 1 : end);
Jm = J(1:m, 1:m);
Jr = J(1:m, m + 1 : end);
dZ = [(Jm \ dJm) * (Jm \ Jr) -  Jm \ dJr; zeros(n - m)];