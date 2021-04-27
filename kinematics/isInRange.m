function [flag, a, b] = isInRange(psi, kesai)
n = length(psi) / 2;
flag = 0;
a = [];
b=[];
for k = 1 : n
    a = psi(2*k-1);
    b = psi(2*k);
    if kesai >= a && kesai <= b
        flag = 1;
        return;
    end
end