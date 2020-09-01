function P_new = ApplyFrameToPoint(f, P, type)
%type == 'angle' or 'rv'
if nargin == 2
    type = 'rv';
end
f = f(:);
r = f(1 : 3);
t = f(4 : 6);
P = reshape(P, 3, []);
N = size(P, 2);
P_new = zeros(3, N);
for i = 1 : N
    P_new(:, i) = ApplyRotationToPoint(r, P(:, i), type) + t;
end