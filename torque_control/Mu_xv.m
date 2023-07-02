function muxv = Mu_xv(J, M, dJ, Z, C)
muxv = (pinv_JT_x(J, M, C) - A_x_x(J, M, dJ)) * Z;