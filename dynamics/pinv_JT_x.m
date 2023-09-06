function pinvJtx = pinv_JT_x(J, M, x)
pinvJtx = (J * (M' \ J')) \ J * (M' \ x);