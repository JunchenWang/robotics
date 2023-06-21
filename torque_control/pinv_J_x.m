function pinvJX = pinv_J_x(J, M, x)
pinvJX = (M \ J') * ((J * (M \ J')) \ x);