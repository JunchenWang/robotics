function muv = Mu_v(Z, M, dZ, dM, C)
muv = (Z' * C - A_v(Z, M) * derivative_pinv_Z(Z, M, dZ, dM)) * Z;