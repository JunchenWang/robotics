function angles =UR_inverse_kin(T, cfg, t6_ref)
%输入关节6相对于基座0的T矩阵
%输出8组可能的解angles，每一列为theta1-6；
%cfg为三个1或-1，分别表示（theta1+fai）的cos值，theta5及theta3的cos正负值，顺序为1,3,5
T60=T;
%UR5
% d1=89.2;
% d4=109.3;
% d5=94.75;
% d6=82.5;
% a3=425;
% a4=392;
% UR5e
d1=162.5;
d4=133.3;
d5=99.7;
d6=99.6;
a3=425;
a4=392.2;
eps0 = 1e-9;%not change
px=d6*T60(1,3)-T60(1,4);
py=T60(2,4)-d6*T60(2,3);
%t1
t1=atan2(-d4/sqrt( px^2+py^2 ),cfg(1)*sqrt(1-(d4^2/( px^2+py^2 )))) -atan2(py,px);
t1 = mod(t1 + pi, 2*pi) - pi;
%t5
t5=cfg(3)*real(acos(T60(1,3)*sin(t1)-T60(2,3)*cos(t1)));
%t6
if abs(t5) <= eps0 || abs(abs(t5)-pi) <= eps0
    t6=t6_ref;
    disp('singularity');
else
    m=T60(2,2)*cos(t1)-T60(1,2)*sin(t1);
    n=T60(2,1)*cos(t1)-T60(1,1)*sin(t1);
    t6=atan2((m*sin(t5))/(m^2+n^2),(-n*sin(t5))/(m^2+n^2));
end
sum234=atan2(-cos(t6)*(T60(1,2)*cos(t1)+T60(2,2)*sin(t1))-sin(t6)*(T60(1,1)*cos(t1)+T60(2,1)*sin(t1)),T60(3,2)*cos(t6)+T60(3,1)*sin(t6));
k1=d5*sin(sum234)-T60(1,4)*cos(t1)+d6*T60(1,3)*cos(t1)+d6*T60(2,3)*sin(t1)-T60(2,4)*sin(t1);
k2=-d5*cos(sum234)-T60(3,4)+d6*T60(3,3)+d1;
%t3
L = k1^2 + k2^2;
if sqrt(L) >= a3 + a4
    error('no solution');
end
t3=cfg(2)*real(acos((L-a4^2-a3^2)/(2*a4*a3)));
%t2
m=a4*cos(t3)+a3;
n=a4*sin(t3);
t2=atan2((k2*m-k1*n)/(m^2+n^2),(k1*m+k2*n)/(m^2+n^2));
%t4
t4=sum234-t3-t2;
t4 = mod(t4 + pi, 2*pi) - pi;
angles=[t1,t2,t3,t4,t5,t6];
end