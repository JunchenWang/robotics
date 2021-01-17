function [loss, J] = trocarAnalysis(v, data)
N = size(data,1);
xyz = v(1:3);
p = v(4:6);
d = v(7:9);
loss = 0;
J = zeros(1, 9);
for i = 1 : N
    Ri = EulZYX2R(data(i, 1:3));
    ti = data(i, 4:6);
    pi = p * Ri' + ti;
    di = d * Ri';
    di2 = di*di';
    mi = dot(xyz-pi,di);
    li2 = (xyz - pi)*(xyz-pi)';
    lossi = li2 - mi*mi / di2;
    loss = loss + lossi;
    Jxyz = [2*(xyz(1) - pi(1))-2*mi*di(1)/di2, 2*(xyz(2) - pi(2))-2*mi*di(2)/di2, 2*(xyz(3) - pi(3))-2*mi*di(3)/di2];
    Jp = -Jxyz*Ri;
    Jd = [-2*mi*(xyz(1)-pi(1))/di2+2*di(1)*mi^2/di2^2, -2*mi*(xyz(2)-pi(2))/di2+2*di(2)*mi^2/di2^2, -2*mi*(xyz(3)-pi(3))/di2+2*di(3)*mi^2/di2^2] * Ri;
    J = J + [Jxyz, Jp, Jd];
end
% penalty = d*d' - 1;
% loss = loss + lambda*(penalty)^2;
% 
% J = J + [0 0 0 0 0 0, 4*lambda*penalty*d];