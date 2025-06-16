function p = BayesianHMMP0(op, x0)
if strcmp(op, 'DS') == 1
    p = 0.3 + sqrt(1) * randn(1, x0);
else
    p = normpdf(x0, 0.3, sqrt(1));
end