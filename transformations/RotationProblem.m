function R = RotationProblem(na, nb)
n = size(na, 2);
na = [zeros(1, n); na];
nb = [zeros(1, n); nb];
A = zeros(4);
for i = 1 : n
    A = A + (Q(na(:, i)) - W(nb(:, i)))' * (Q(na(:, i)) - W(nb(:, i)));
end
[V D] = eig(A);
R = RMByQuaternion(V(:, 1));