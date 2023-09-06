function Axx = A_x_x(J, M, x)
Axx = (J * (M \ J')) \ x;