function p = BayesianHMM(op, xk, xk_1, w)
if strcmp(op, 'DS') == 1
    p = xk_1 + sqrt(w) * randn(1, xk);
else
    p = normpdf(xk - xk_1, 0, sqrt(w));
end