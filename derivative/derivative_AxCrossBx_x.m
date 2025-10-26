function J = derivative_AxCrossBx_x(A, B, x)
a1 = A(:, 1);
a2 = A(:, 2);
a3 = A(:, 3);

b1 = B(:, 1);
b2 = B(:, 2);
b3 = B(:, 3);

C = [cross(a1, b1), cross(a1, b2), cross(a1, b3),...
     cross(a2, b1), cross(a2, b2), cross(a2, b3),...
     cross(a3, b1), cross(a3, b2), cross(a3, b3)];
D = [2 * x(1) 0 0
     x(2) x(1) 0
     x(3) 0 x(1)
     x(2) x(1) 0
     0 2 * x(2) 0
     0 x(3) x(2)
     x(3) 0 x(1)
     0 x(3) x(2)
     0 0 2 * x(3)];
 J = C * D;
     