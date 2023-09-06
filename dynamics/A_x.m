function Ax = A_x(J, M)
Ax = inv(J * (M \ J'));