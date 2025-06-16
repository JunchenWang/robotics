function p = BayesianLikelihood(yk, xk, v)
p = normpdf(yk - xk, 0, sqrt(v));
