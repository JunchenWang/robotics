function dist = twist_dist(T1, T2)
% 计算T1到T2两个SE矩阵的距离 T1*T(dist) = T2
dist = logT(InvertT(T1)*T2);