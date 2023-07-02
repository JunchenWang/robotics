function muxx = Mu_x_x(J, M, dJ, C, x)
muxx = (pinv_JT_x(J, M, C) - A_x_x(J, M, dJ)) * pinv_J_x(J, M, x);