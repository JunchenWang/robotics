function [mu_x_t sigma_x_t mu_y_t sigma_y_t sigma_xy_t] = KalmanPrediction(ft, ht, Rt, Qt, mu_x, sigma_x, Cx_t, Cy_t, type, Params)
if nargin < 9
    type = 'L';
end
if nargin < 10
    Params = [];
end

if strcmp(type, 'L')
    if isempty(Cx_t)
        mu_x_t = ft * mu_x;
    else
        mu_x_t = ft * mu_x + Cx_t;
    end
    if isempty(Cy_t)
        mu_y_t = ht * mu_x_t;
    else
        mu_y_t = ht * mu_x_t + Cy_t;
    end
    sigma_x_t = ft * sigma_x * ft' + Qt;
    sigma_y_t = ht * sigma_x_t * ht' + Rt;
    sigma_xy_t = sigma_x_t * ht';
    return;
end

if strcmp(type, 'E')
    Dim_R = size(Rt, 1);
    Dim_Q = size(Qt, 1);
    if isempty(Cx_t)
        [mu_x_t A W] = ft(mu_x, zeros(Dim_Q, 1)); %#ok<RHSFN>
    else
        [mu_x_t A W] = ft(mu_x, zeros(Dim_Q, 1), Cx_t);%#ok<RHSFN>
    end
    if isempty(Cy_t)
        [mu_y_t H V] = ht(mu_x_t, zeros(Dim_R, 1));%#ok<RHSFN>
    else
        [mu_y_t H V] = ht(mu_x_t, zeros(Dim_R, 1), Cy_t);%#ok<RHSFN>
    end
    sigma_x_t = A * sigma_x * A' + W * Qt * W';
    sigma_y_t = H * sigma_x_t * H' + V * Rt * V';
    sigma_xy_t = sigma_x_t * H';
    return;
end

if strcmp(type, 'U')
    Dim_R = size(Rt, 1);
    Dim_Q = size(Qt, 1);
    N = length(mu_x);
    if isempty(Cx_t)
        ft = @(x) ft(x(1 : N), x(N + 1 : N + Dim_Q));
    else
        ft = @(x) ft(x(1 : N), x(N + 1 : N + Dim_Q), Cx_t);
    end
    if isempty(Cy_t)
        ht = @(x) ht(ft(x), x(N + Dim_Q + 1 : N + Dim_Q + Dim_R)); %要注意ht(ft(x), ...)
    else
        ht = @(x) ht(ft(x), x(N + Dim_Q + 1 : N + Dim_Q + Dim_R), Cy_t);%要注意ht(ft(x), ...)
    end
    UT_mu = [mu_x; zeros(Dim_Q, 1); zeros(Dim_R, 1)];
    UT_sigma = [sigma_x, zeros(N, Dim_Q), zeros(N, Dim_R);
                zeros(Dim_Q, N), Qt, zeros(Dim_Q, Dim_R); 
                zeros(Dim_R, N), zeros(Dim_R, Dim_Q), Rt];
     if isempty(Params)
         [mu_x_t sigma_x_t mu_y_t sigma_y_t sigma_xy_t] = UT_Cross(UT_mu, UT_sigma, ft, ht);
     elseif length(Params) == 1
         [mu_x_t sigma_x_t mu_y_t sigma_y_t sigma_xy_t] = UT_Cross(UT_mu, UT_sigma, ft, ht, Params(1));
     elseif length(Params) == 2
         [mu_x_t sigma_x_t mu_y_t sigma_y_t sigma_xy_t] = UT_Cross(UT_mu, UT_sigma, ft, ht, Params(1), Params(2));
     else
         [mu_x_t sigma_x_t mu_y_t sigma_y_t sigma_xy_t] = UT_Cross(UT_mu, UT_sigma, ft, ht, Params(1), Params(2), Params(3));
     end
     return;
end

error('Undefined KalmanPrediction Method!');