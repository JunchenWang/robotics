function test_LMI
% 李雅普诺夫方程
% ATX + XA < 0, -X < 0
n = 6;
A = rand(n) - 10*eye(n);
setlmis([]);
X = lmivar(1,[n, 1]);
lmiterm([1 1 1 X], 1, A, 's'); %XA + ATX
lmiterm([2 1 1 X], -1, 1); % -X
lmisys = getlmis;
[tmin,xfeas] = feasp(lmisys);
X = dec2mat(lmisys, xfeas, X);
disp(X);
eig(X)

