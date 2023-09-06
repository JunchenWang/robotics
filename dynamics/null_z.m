function Z = null_z(J)
% each col of Z is the null vector of J
% J has full row rank
[m, n] = size(J);
Jm = J(1:m, 1:m);
Jr = J(1:m, m + 1 : end);
Z = [-Jm \ Jr; eye(n - m)];
