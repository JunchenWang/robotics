function muvx = Mu_vx(J, M, dZ, dM, Z, C)
muvx = (Z' * C - A_v(Z, M) * derivative_pinv_Z(Z, M, dZ, dM)) * pinv_J(J, M);