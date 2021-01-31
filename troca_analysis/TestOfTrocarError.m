%一排包括两个小球的相机坐标， 距离法兰远的小球排在前面。
%第一排是开始前采集，使用第一个点作为trocar点
p =[
21.763 -81.091 -1288.659 22.497 1.268 -1319.648
21.255 -81.717 -1288.273 40.197 -2.664 -1322.014
22.494 -81.602 -1288.705 -13.892 -5.782 -1314.550
21.573 -82.805 -1289.938 30.254 3.483 -1275.169
22.530 -84.343 -1284.940 10.519 -20.751 -1344.855
18.547 -191.393 -1275.545 21.440 -104.186 -1286.946



];


B =p(:,4:6); 
A= p(:,1:3);
%第一排是开始前采集，使用第一个点作为trocar点
trocarP=[p(1,1), p(1,2),p(1,3)];

errors = zeros (1,size(p,1));
for n = 1:size(p,1)

    
errors(n) = trocarPointDistanceToLineOf2Points(trocarP, A(n,:), B(n,:));

end


%trocar 点到电切镜轴的距离
function d = trocarPointDistanceToLineOf2Points(trocarPoint, pointA,pointB)
    ta = trocarPoint -pointA;
    ba = pointB - pointA;
    d =norm( cross(ta, ba))/norm(ba);
    
   
end