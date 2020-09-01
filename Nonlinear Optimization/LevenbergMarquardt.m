function [xk Residual g_norm] = LevenbergMarquardt(f, x0, epsilon, XTOL, FTOL, delta_k, state)
if nargin < 3
    epsilon = [];%1e-16;
end
if nargin < 4
    XTOL = [];%1e-16;
end
if nargin < 5
    FTOL = [];%1e-16;
end
if nargin < 6
    delta_k = [];%norm(x0);
end
if nargin < 7
    state = 1;
end
if isempty(epsilon)
    epsilon = 1e-16;
end
if isempty(XTOL)
    XTOL = 1e-16;
end
if isempty(FTOL)
    FTOL = 1e-16;
end
if isempty(delta_k)
    delta_k = norm(x0);
end
 
f_limit = 5000;
iter_limit = 5000;
xk = x0;
[fk Jk] = f(xk);
f_count = 1;
iter_n = 0;
sigma = 0.1;
[m n] = size(Jk);
Dk = zeros(n);
for i = 1 : n
    Dk(i, i) = norm(Jk(:, i));
end
g_norm = norm(Jk' * fk);
f_norm = norm(fk);
if state
    display('****************************** Begin Levenberg-Marquardt ******************************');
    msg = sprintf('epsilon: %0.1e\t XTOL: %0.1e\t FTOL: %0.1e\t f_limit: %d\t iter_limit: %d\t delta_k: %0.5e\t g_norm: %0.5e\t f_norm: %0.5e',...
        epsilon, XTOL, FTOL, f_limit, iter_limit, delta_k, g_norm, f_norm);
    disp(msg);
    disp(' ');
    finish_str = {'epsilon', 'XTOL', 'f_count', 'iter_n', 'FTOL'};
    t_begin = tic;
end

while 1
    
    if g_norm <= epsilon
        flag = 1;
        break;
    end
    if delta_k <= norm(Dk * xk) * XTOL
        flag = 2;
        break;
    end
    if f_count > f_limit
        flag = 3;
        break;
    end
    if iter_n > iter_limit
        flag = 4;
        break;
    end
    
    pk = -pinv(Jk) * fk;
    if norm(Dk * pk) <= (1 + sigma) * delta_k
        lambda_k = 0;
    else
        lambda_k = get_alpha(fk, Jk, Dk, delta_k, sigma);
        pk = determine_p(fk, Jk, Dk, lambda_k);
    end
    [rho_k fn ffn tem1 tem2] = get_rho(xk, fk, Jk, pk, Dk, f, lambda_k);
    f_count = f_count + 1;
    if rho_k > 0.0001
        xk = xk + pk;
        iter_n = iter_n + 1;
        [fk Jk] = f(xk);
        g_norm = norm(Jk' * fk);
        f_norm = norm(fk);
        step_norm = norm(pk);
        f_count = f_count + 1;
        if state
            msg = sprintf('iter_n: %3d\t f_count: %3d\t f_norm: %0.5e\t g_norm: %0.5e\t step_norm: %0.5e\t delta_k: %0.5e\t lambda_k: %0.5e',...
                iter_n, f_count, f_norm, g_norm, step_norm, delta_k, lambda_k);
            disp(msg);
        end
    end
    
    if tem1 + 2 * lambda_k * tem2 <= FTOL
        flag = 5;
        break;
    end
    
    if rho_k <= 0.25
        mu = get_mu(fn, ffn, tem1, tem2, lambda_k);
        delta_k = mu * delta_k;
        if state
            msg = sprintf('reduce delta_k: %0.5e', delta_k);
            disp(msg);
        end
    elseif (rho_k > 0.25 && rho_k < 0.75 && lambda_k == 0) ...
           || (rho_k >= 0.75)
        delta_k = 2 * norm(Dk * pk);
        if state
            msg = sprintf('enlarge delta_k: %0.5e', delta_k);
            disp(msg);
        end
    end
    for i = 1 : n
        Dk(i, i) = max(norm(Jk(:, i)), Dk(i, i));
    end
end
Residual = sqrt(f_norm^2 / m);
if state
    disp(' ');
    toc(t_begin);
    msg = sprintf('Conditions: %s\t Residual: %f', finish_str{flag}, Residual);
    disp(msg);
    display('****************************** Levenberg-Marquardt End ******************************');
end

function [pk Ek R_lambda] = determine_p(fk, Jk, Dk, lambda_k)
[m n] = size(Jk);
[Qk Rk Ek] = qr(Jk);
Qk = Qk';
D_lambda = sqrt(lambda_k) * Ek' * Dk * Ek;
[W RR] = qr([Rk; D_lambda], 0);
W = W';
R_lambda = RR(:, 1 : n);
tem = zeros(m + n, 1);
tem(1 : m) = Qk * fk;
u = W * tem;
pk = -Ek / R_lambda * u;


function [rho fn ffn tem1 tem2] = get_rho(xk, fk, Jk, pk, Dk, f, lambda_k)
ffn = norm(f(xk + pk));
fn = norm(fk);
Jpn = norm(Jk * pk);
Dpn = norm(Dk * pk);
tem1 = (Jpn / fn)^2;
tem2 = (Dpn / fn)^2;
if ffn > fn
    rho = 0;
else
    rho = (1 - (ffn / fn)^2) / (tem1 + 2 * lambda_k * tem2);
end

function mu = get_mu(fn, ffn, tem1, tem2, lambda_k)
if ffn <= fn
    mu = 0.5;
elseif ffn <= 10 * fn
    gamma = -(tem1 + lambda_k * tem2);
    mu = gamma / (2 * gamma + 1 - (ffn / fn)^2);
else
    mu = 0.1;
end

function alpha = get_alpha(fk, Jk, Dk, delta_k, sigma)
[m n] = size(Jk);
JJ = Jk / Dk;
[U D V] = svd(JJ, 0); %#ok<NASGU>
z = U' * fk;
s = diag(D);
tol = max(m, n) * eps(max(s));
r = sum(s > tol);
phi = @(alpha) norm(s .* z(1 : n) ./ (s.^2 + alpha)) - delta_k;
phi_d = @(alpha) phi_derive(alpha, fk, Jk, Dk);
u0 = norm((Jk / Dk)' * fk) / delta_k;
if r == n
    l0 = -phi(0) / phi_d(0);
else
    l0 = 0;
end
alpha = u0;
while 1%abs(phi(alpha)) > sigma * delta_k
%     display('trying to determine alpha');
    if alpha < l0 || alpha > u0
        alpha = max(0.001 * u0, sqrt(l0 * u0));
    end
    if abs(phi(alpha)) <= sigma * delta_k
        break;
    end
    phi_alpha = phi(alpha);
    phi_d_alpha = phi_d(alpha);
    if phi_alpha < 0
        u0 = alpha;
    end
    l0 = max(l0, alpha - phi_alpha / phi_d_alpha);
    alpha = alpha - (phi_alpha + delta_k) / delta_k * phi_alpha / phi_d_alpha;
end


function phi_d = phi_derive(alpha, fk, Jk, Dk)
[p_alpha E R_alpha] = determine_p(fk, Jk, Dk, alpha);
q = Dk * p_alpha;
qn = norm(q);
phi_d = -qn * (norm(R_alpha' \ (E' * Dk' * q) / qn))^2;
