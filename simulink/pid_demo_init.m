Ts = 0.01;
Tf = 0.02;
kp = 0.5;
ki = 2;
kd = 10;
A = [0, 0, 0; 
    -kd / (Tf + Ts), Tf / (Tf + Ts), 0; 
    0, 0, 1];
B = [1; kd / (Tf + Ts); ki * Ts];
C = [-kd / (Tf + Ts), Tf / (Tf + Ts), 1];
D = kp + kd / (Tf + Ts);

