function dh_row = rttr_inverse(T)
% dh_row: alpha, a, d, theta
% T: 4x4 Matrix,��������T���ж�Ӧ��DH����
theta = atan2(-T(1,2), T(1,1));
alpha = atan2(-T(2,3), T(3,3)); 
a = T(1,4);
spha = sin(alpha);
cpha = cos(alpha);
if abs(spha) > abs(cpha)
    d = T(2,4) / -spha;
else
    d = T(3,4) / cpha;
end
dh_row = [alpha, a, d, theta];

