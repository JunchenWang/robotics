function t = test_LMI_DO
% 李雅普诺夫方程
% ATX + XA < 0, -X < 0
n = 7;
m = 20;
zetas = [1.5, 2, 2.5, 3, 3.5, 4];
K = length(zetas);

t = zeros(K,m);
y = linspace(2, 10, m);

for k = 1 : K
    zeta = zetas(k);
for i = 1 : m
Y = y(i) * eye(7);

setlmis([]);
% Y = lmivar(1,[n, 0]);
Gamma = lmivar(1,[n, 1]);

lmiterm([-1 1 1 0], Y +Y'-zeta*eye(n)); %Y + YT
% lmiterm([-1 1 1 Y], 1, 1, 's'); %Y + YT
lmiterm([-1 1 1 0], -zeta*eye(n)); %-zeta I
lmiterm([-1 1 2 0], Y'); %YT
% lmiterm([-1 1 2 -Y], 1,1); %YT
lmiterm([-1 2 2 Gamma], 1,1); %Gamma

% lmiterm([-2 1 1 Y], 1, 1); % 0 < Y
lmiterm([-2 1 1 Gamma], 1, 1); % 0 < Gamma
lmisys = getlmis;
[tmin,xfeas] = feasp(lmisys);
t(k,i) = tmin;
end
end
plot(y, t, 'LineWidth',2);
title('$\zeta=2$','interpreter','latex');
xlabel("$y$", 'interpreter','latex');
ylabel('best $\sigma$', 'interpreter','latex','Rotation',0);
yticks([-20, -15, -10, -5, 0, 5]);
xticks([2,3,4,5,6,7,8,9,10]);
set(gca,'FontSize', 32);
set(gcf,'Position',[100 100 1200 800]);
grid on;
hold on;
plot(y, zeros(1, length(y)), 'black-.', 'LineWidth',2);
hold off;
legend('$\zeta=1.5$','$\zeta=2.0$','$\zeta=2.5$','$\zeta=3.0$','$\zeta=3.5$','$\zeta=4.0$','interpreter','latex');
% fontsize(lg,18,'points')
% Y = dec2mat(lmisys, xfeas, Y);
% Gamma = dec2mat(lmisys, xfeas, Gamma);
% disp(Y);
% disp(Gamma);

