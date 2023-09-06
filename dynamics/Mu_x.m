function mux = Mu_x(J, M, dJ, C)
mux = (pinv_JT_x(J, M, C) - A_x_x(J, M, dJ)) * pinv_J(J, M);