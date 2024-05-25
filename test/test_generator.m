fid = fopen('E:\manuscripts\孙振\完整RCM数据&计算结果\p.txt','r');
data = fscanf(fid, '%f %f %f %f', [4 inf]);
data(1,:) = data(1,:) / 1000;
color='rgbcykm';
K = size(data,2);
t = 1 : K;
t = t * 2e-3;
t(t>50) = [];
K = length(t);
for i = 1 : 4
    plot(t, data(i,1:K), color(mod(i - 1,7)+1));
    hold on;
end
fclose(fid);
figure
fid = fopen('E:\manuscripts\孙振\完整RCM数据&计算结果\trajectory.txt','r');
data = fscanf(fid, '%f %f %f %f %f %f %f', [7 inf]);
for i = 1 : 7
    plot(t, data(i,1:K), color(mod(i - 1,7)+1));
    hold on;
end
fclose(fid);

figure;
fid = fopen('E:\manuscripts\孙振\完整RCM数据&计算结果\omegadata.txt','r');
data = fscanf(fid, '%f %f %f %f', [4 inf]);
data(1,:) = data(1,:) / 1000;
K = size(data,2);
t = 1 : K;
t = t * 1e-2;
t(t>50) = [];
K = length(t);
for i = 1 : 4
    plot(t, data(i,1:K), color(mod(i - 1,7)+1));
    hold on;
end
fclose(fid);
% return;
% u = udpport("byte");
% Freq = 250;
% r = rateControl(Freq);
% for i = 1 : K
%     cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', data(1,i), data(2,i), data(3,i), data(4,i), data(5,i)...
%     ,data(6,i), data(7,i));
%     writeline(u,cmd,"192.168.3.34",7755);
%     waitfor(r);
% end
