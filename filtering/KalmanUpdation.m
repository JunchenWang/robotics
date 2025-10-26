function  [mu_x sigma_x] = KalmanUpdation(mu_x_t, sigma_x_t, mu_y_t, sigma_y_t, sigma_xy_t, yt)
 K = sigma_xy_t / sigma_y_t;
 y_err_t = yt - mu_y_t;
 mu_x = mu_x_t + K * (y_err_t);
 sigma_x = sigma_x_t - K * sigma_y_t * K';