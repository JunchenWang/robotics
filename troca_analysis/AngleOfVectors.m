euler1=[0 0 0];

arr = zeros(37,2);

for i=10:10:360 
euler2=[i 0 0];


euler1 = euler1/180 * pi;
euler2 = euler2/180 * pi;



Q1=angle2quat(euler1(1),euler1(2),euler1(3),  'ZYX');

Q2=angle2quat(euler2(1),euler2(2),euler2(3),  'ZYX');


conjQ1 = quatconj(Q1);
Q12 = quatmultiply(conjQ1,Q2);

angle = 2 * atan2(norm(Q12(2:4)),Q12(1)) /pi*180;

angleB = 2* acos(Q12(1))/pi*180;
arr( floor(i/10) ,:) = [Q12(1) angleB];
end