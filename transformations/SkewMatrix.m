function M = SkewMatrix(x)
if isvector(x)
    M = [0 -x(3) x(2);
        x(3) 0  -x(1);
        -x(2) x(1) 0];
else
    M = [x(3, 2); x(1, 3); x(2, 1)];
end
 