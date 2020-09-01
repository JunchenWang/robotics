function qd2 = qd_converter(qd1)
if isvector(qd1)
    qd2 = reshape(qd1, 4, 2);
else
    qd2 = qd1(:);
end