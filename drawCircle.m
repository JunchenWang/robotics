function [] = drawCircle()
tipRangeOfScope();
t = 0:pi/500:pi;
ir = [
    0 0
    10 4
    20 8
    30 12
    40 16
    50 20
    60 24
    70 28
    80 32
    90 36
    100 40
    110 42.5
    120 45
    130 47.5
    140 50  
    150 52.5
    160 55
    170 57.5
    180 60
    190 62.5
    200 64
    210 65
    220 60
    230 48
    240 0

    
    ];
drift = [-995.7579756435572 10.862778948329678 39.91170326096567];

for i =1:size(ir,1)
    zt(1,:) = sin(2*t)*ir(i,2)  + drift(3);
    yt(1,:) = cos(2*t) *ir(i,2) + drift(2);
    xt = t-t + drift(1) - ir(i,1) ;


    plot3(xt,yt,zt,'ro');
    hold on;
end 


   
%  [x,y,z]=cylinder([40,40],40);
%   surf( -180 * z  + drift(1)-100,y+ drift(2),x+ drift(3));
  
%  [x,y,z]=cylinder([50,50],40);
%  surf( -130 * z  + drift(1)-140,y+ drift(2),x+ drift(3));
%  
%   [x,y,z]=cylinder([0,50],40);
%  surf( -140 * z  + drift(1),y+ drift(2),x+ drift(3));
 
  [x,y,z]=cylinder(ir(:,2)',40);
   surf( -240 * z  + drift(1),y+ drift(2),x+ drift(3));

axis equal;

irWithangle = zeros (size(ir,1),3);
irWithangle(:,1:2) = ir;
for n = 2: size(ir,1)

    angle = ir(n,2)/ir(n,1);
    irWithangle(n,3) =    atan(angle)/pi*180;
    
end
irWithangle
end



