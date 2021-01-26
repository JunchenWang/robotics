%body frames
T =[
    15.9352307 -78.9342745 -82.2367877 -31.749 206.835 -1359.300
42.8779307 -76.3448401 -145.1530702 40.693 214.680 -1385.914
-69.3049201 -67.2794789 -28.3510074 35.619 201.216 -1321.858
86.2239527 -63.4523132 -172.3622602 5.171 211.121 -1416.690
-76.3254372 -88.5797215 -18.9398131 16.725 114.263 -1361.965


];

% coordinates on body frame.
p=[
-21.291 -11.042 116.669
-21.031 9.685 114.475
-21.739 -0.784 -123.614

];
B = (p(1,:) + p(2,:) )/2;
A= p(3,:);
% camera coordinates
trocarP=[12.726	92.643	-1383.735]';

errors = zeros (1,size(T,1));
for n = 1:size(T,1)
    Trans =rot2t( T(n,:));
    Bcamera = Trans * [B 1]';
    Acamera = Trans * [A 1]';
    
errors(n) = trocarPointDistanceToLineOf2Points(trocarP, Acamera(1:3,1), Bcamera(1:3,1));

end

% NDI数据  欧拉角转为旋转矩阵加上平移
function T = rot2t(data)
    
zyx = data(1:3);
zyx = zyx/180*pi;
t = data(4:6);
t= t';
R = eul2rotm(zyx,'zyx');
T = [R,t;0 0 0 1];
end

%trocar 点到电切镜轴的距离
function d = trocarPointDistanceToLineOf2Points(trocarPoint, pointA,pointB)
    ta = trocarPoint -pointA;
    ba = pointB - pointA;
    d =norm( cross(ta, ba))/norm(ba);
    
   
end