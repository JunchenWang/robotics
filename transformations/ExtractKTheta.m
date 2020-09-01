function [k, theta] = ExtractKTheta(f)
r = f(1 : 3);
r =Normalize_r(r);
theta = norm(r);
k = r / theta;
