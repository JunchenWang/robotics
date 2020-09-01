function qdn = norm_qd(qd)
qd0n = norm(qd(1 : 4));
qdn = [qd0n, dot(qd(1 : 4), qd(5 : 8)) / qd0n]'; 