function pole = cal_second_order_pole(Mp, Ts)
% 对于欠阻尼系统给出超调量Mp和调节时间Ts，计算共轭极点
a = log(Mp);
zeta = sqrt(a.^2 / (pi.^2 + a.^2));
omega_n = 4 / (Ts * zeta);
pole = [-zeta * omega_n + omega_n * sqrt(1 - zeta.^2) * 1i, -zeta * omega_n - omega_n * sqrt(1 - zeta.^2) * 1i];