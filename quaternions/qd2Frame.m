function f = qd2Frame(qd)
n = size(qd, 2);
f = zeros(6, n);
for i = 1 : n
    q = qd(1:4, i);
    qe = qd(5:8, i);
    r = r_q_converter(q);
    t = 2 * ConcatenateQuaternions(qe, Conjugate_q(q));
    f(:, i) = [r;t(2:4)];
end