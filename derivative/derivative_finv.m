function J = derivative_finv(f)
f = f(:);
r = f(1 : 3);
t = f(4 : 6);
J = [-eye(3) zeros(3)
      derivative_rx_r(-r, t) -RotationByAxisAngleRep(r)'];