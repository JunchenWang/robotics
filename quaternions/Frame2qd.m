function qd = Frame2qd(f)
n = size(f, 2);
qd = zeros(8, n);

for i = 1 : n
    r = f(1:3, i);
    t = f(4:6, i);
    q = r_q_converter(r);
    qd(:, i) = [q; 0.5 * ConcatenateQuaternions([0; t], q)];
end