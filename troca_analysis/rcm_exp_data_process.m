fid = fopen("E:\manuscripts\RCM-control\data\240503_afternoon\240503_afternoon\fast_rcmvf_Y2_1.txt");
data = fscanf(fid, '%f', [29, inf]);
time_s = 10;%起始时间
time_e = 50;%终止时间
t = data(1,:);
id = (t >= time_s) & (t <= time_e);
t = t(id) - time_s;
xe = data(2:4,id);
dq = data(9:15,id);
tau_d = data(16:22,id);
tau_c = data(23:29,id);

fs = 24;
ls = 24;
tl = tiledlayout(2,2);
nexttile
plot(t,sqrt(sum(xe.^2,1)) * 1000, 'LineWidth',2);
grid on;
% xlabel("$t$/s", 'interpreter','latex');
ylabel('$||x_e||$/mm', 'interpreter','latex');
set(gca,'FontSize', fs);
yticks([0, 0.5, 1, 1.5]);
% xticks([0,1,2,3,4,5,6,7,8,9,10]);

nexttile
plot(t,dq,'LineWidth',2);
grid on;
% xlabel("$t$/s", 'interpreter','latex');
ylabel('$\dot{q}$/(rad/s)', 'interpreter','latex');
set(gca,'FontSize', fs);
lg = legend('J1', 'J2','J3','J4','J5','J6','J7', 'Orientation','horizontal');
fontsize(lg,ls,'points')


nexttile
plot(t,tau_d,'LineWidth',2);
grid on;
% xlabel("$t$/s", 'interpreter','latex');
ylabel('$\hat{\tau}_x$/(Nm)', 'interpreter','latex');
set(gca,'FontSize', fs);
lg = legend('J1', 'J2','J3','J4','J5','J6','J7','Orientation','horizontal');
fontsize(lg,ls,'points')

nexttile
plot(t,tau_c,'LineWidth',2);
grid on;
% xlabel("$t$/s", 'interpreter','latex');
ylabel('$\tau$/(Nm)', 'interpreter','latex');
set(gca,'FontSize', fs);
lg = legend('J1', 'J2','J3','J4','J5','J6','J7', 'Orientation','horizontal');
fontsize(lg,ls,'points')

tl.TileSpacing = 'tight';
tl.Padding = 'compact';
title(tl,'Fast Interaction with NDOB', 'FontSize', fs);
xlabel(tl,"$t$/s", 'interpreter','latex', 'FontSize', fs);
fclose(fid);