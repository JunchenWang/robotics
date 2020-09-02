function tf = RVFrame2EulZYXFrame(rvf)
n = size(rvf, 2);
tf = rvf;
for i = 1 : n
    tf(1:3, i) = R2EulZYX(exp_w(rvf(1:3,i)));
end