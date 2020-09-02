function rvf = EulZYXFrame2RVFrame(tf)
n = size(tf, 2);
rvf = tf;
for i = 1 : n
    rvf(1:3, i) = logR(EulZYX2R(tf(1:3,i)));
end
