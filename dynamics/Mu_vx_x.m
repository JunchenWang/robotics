function muvxx = Mu_vx_x(J, Z, M, dZ, dM, C, x)
muvxx = (Z' * C - A_v(Z, M) * derivative_pinv_Z(Z, M, dZ, dM)) * pinv_J_x(J, M, x);