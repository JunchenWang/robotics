%相机坐标系
    %八个螺丝坐标，
    %工具小球A
    %工具小球C
p =[

-39.676 313.348 -1371.709
-11.092 320.779 -1352.623
23.504 321.384 -1358.409
44.041 314.429 -1386.270
38.520 304.224 -1419.443
10.554 296.724 -1438.567
-24.409 296.361 -1432.945
-44.952 303.027 -1405.280

21.361	-79.899	-1289.595	%OK	48.186	-40.903	-1306.334	OK

18.046	4.391	-1314.579	%OK	-23.99	-40.344	-1293.491

];

% len的长度应该大概相同。 用以校验NDI采集螺丝坐标数据的准确性
len = zeros(1,7);
for n =1:7
    len (n) = norm (p(n,:) - p(n+1,:));
end


p1= (p(8,:) + p(1,:) )/2;
p2= (p(2,:) + p(3,:) )/2;
p3= (p(4,:) + p(5,:) )/2;
p4= (p(6,:) + p(7,:) )/2;
o =(p1 + p2 + p3+ p4)/4;


x= p1 -p3;
y=p4-p2;
z= cross (x,y);
y = cross (z,x);

x = x/norm(x);
y = y/norm(y);
z = z/norm(z);
R = [x',y',z'];


T = [R,o';0 0 0 1];
R = eye(3);
o = [0 0 24];
T = T * [R,o';0 0 0 1];


validateP = zeros(size(p,1),4 );
for n = 1: size(p,1)
validateP(n,:) =  (T\    [ p(n,:)';1])';
end

validateP = validateP(:,1:3);
trocarPoint = validateP(9,:);
direction = validateP(9,:) - validateP(10,:)     ;
direction = direction / norm(direction);

v= validateP(1:8,:);
figure
plot3(v(:,1),v(:,2),v(:,3))



