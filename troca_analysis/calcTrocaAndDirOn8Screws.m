%八个螺丝坐标，三个电切镜轴坐标。 最后一个同时也是不动点坐标
p =[
    
-53.670 306.571 -1594.568
-25.456 311.146 -1574.375
9.245 311.075 -1579.675
30.247 306.996 -1607.687
25.021 300.715 -1642.021
-3.162 296.029 -1662.564
-38.139 296.233 -1656.981
-58.979 300.444 -1628.687
-19.238 15.384 -1595.764
-0.331 14.234 -1598.277
-7.163 -219.709 -1554.415


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
trocarPoint = validateP(11,:);
direction = validateP(11,:) - (validateP(9,:) + validateP(10,:) )/2    ;
direction = direction / norm(direction);

v= validateP(1:8,:);
figure
plot3(v(:,1),v(:,2),v(:,3))



