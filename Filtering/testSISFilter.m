function testSISFilter
T = 100;
R = 0.2^2;
Q = 0.01^2;
P0 = 0.1^2;
x0 = 1;
xt = x0 + sqrt(Q) * randn(1, T);
y = xt + sqrt(R) * randn(1, T);
x2 = KalmanFilter(1, 1, x0, P0, y, R, Q);
x = SISBayesianFilter(@prior, @likelihood, x0, P0, y, R, Q, @proposal);
plot([x0 xt],'b');
hold on;
plot(x2,'g');
plot(x,'c');
plot(2 : T + 1, y, 'r');

function p = prior(op, xk1, xk, Q)
if strcmp(op, 'p')
    p = normpdf(xk1, xk, sqrt(Q));
else
    p = xk1 + xk;
end
function p = likelihood(op, yk, xk, R)
if strcmp(op, 'p')
    p = normpdf(yk, xk, sqrt(R));
else
    p = yk + xk;
end

function [sample p Pt] = proposal(ft, ht, xt_1, Pt_1, yt, Rt, Qt, Cx_t, Cy_t)
    Dim_R = size(Rt, 1);
    Dim_Q = size(Qt, 1);
    N = length(xt_1);
    if isempty(Cx_t)
        ft = @(x) ft('e', x(1 : N), x(N + 1 : N + Dim_Q));
    else
        ft = @(x) ft('e', x(1 : N), x(N + 1 : N + Dim_Q), Cx_t);
    end
    if isempty(Cy_t)
        ht = @(x) ht('e', ft(x), x(N + Dim_Q + 1 : N + Dim_Q + Dim_R)); %要注意ht(ft(x), ...)
    else
        ht = @(x) ht('e', ft(x), x(N + Dim_Q + 1 : N + Dim_Q + Dim_R), Cy_t);%要注意ht(ft(x), ...)
    end
    UT_mu = [xt_1; zeros(Dim_Q, 1); zeros(Dim_R, 1)];
    UT_sigma = [Pt_1, zeros(N, Dim_Q), zeros(N, Dim_R);
                zeros(Dim_Q, N), Qt, zeros(Dim_Q, Dim_R); 
                zeros(Dim_R, N), zeros(Dim_R, Dim_Q), Rt];
   
    [mu_x_t sigma_x_t mu_y_t sigma_y_t sigma_xy_t] = UT_Cross(UT_mu, UT_sigma, ft, ht, 1 , 0 , 2);
    [mu_x sigma_x] = KalmanUpdation(mu_x_t, sigma_x_t, mu_y_t, sigma_y_t, sigma_xy_t, yt);
    sample = mvnrnd(mu_x', sigma_x)';
    p = mvnpdf(sample', mu_x', sigma_x);
    Pt = sigma_x;


function xt = f(xt_1,w)
xt = xt_1 + w;

function yt = h(xt,v)
yt = xt + v;