%%  plant model
A = -diag([1, 2]);
B = [3; 4];
C = eye(2);
D = 0;
plant = ss(A, B, C, D);
step(plant);




