s = tf('s');
zeta = -1;
omega_n = 10;
G = omega_n.^2 / (s^2 + 2 * zeta * omega_n * s + omega_n.^2);
step(G, [0, 1]);