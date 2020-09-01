function [xk Residual g_norm] = Gauss_Newton(f, x0, epsilon)
if nargin < 3
    epsilon = 1e-16;
end

f_limit = 5000;
iter_limit = 5000;
xk = x0;
[fk Jk] = f(xk);
f_count = 1;
iter_n = 0;
m = size(Jk, 1);
g_norm = norm(Jk' * fk);
f_norm = norm(fk);
display('****************************** Begin G-N ******************************');
msg = sprintf('epsilon: %0.1e\t f_limit: %d\t iter_limit: %d\t g_norm: %0.5e\t f_norm: %0.5e',...
              epsilon, f_limit, iter_limit, g_norm, f_norm);
disp(msg);
disp(' ');
finish_str = {'epsilon', 'No_descent', 'f_count', 'iter_n'};
t_begin = tic;
while 1
    
    if g_norm <= epsilon
        flag = 1;
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
    ffn = norm(f(xk + pk));
    fn = norm(fk);
    Jpn = norm(Jk * pk);
    tem1 = (Jpn / fn)^2;
    if ffn > fn
        rho_k = 0;
    else
        rho_k = (1 - (ffn / fn)^2) / tem1;
    end
    f_count = f_count + 1;
    if rho_k > 0
        xk = xk + pk;
        iter_n = iter_n + 1;
        [fk Jk] = f(xk);
        g_norm = norm(Jk' * fk);
        f_norm = norm(fk);
        step_norm = norm(pk);
        f_count = f_count + 1;
        msg = sprintf('iter_n: %3d\t f_count: %3d\t f_norm: %0.5e\t g_norm: %0.5e\t step_norm: %0.5e',...
            iter_n, f_count, f_norm, g_norm, step_norm);
        disp(msg);
    else
        flag = 2;
        break;
    end
end
disp(' ');
toc(t_begin);
Residual = sqrt(f_norm^2 / m);
msg = sprintf('Conditions: %s\t Residual: %f', finish_str{flag}, Residual);
disp(msg);
display('****************************** G-N End ******************************');
