function x = SISBayesianFilter(f, h, x0, P0, y, R, Q, ImportanceFunc, Ns, IsDebug, Gt, Cx, Cy)
error(nargchk(8, 12, nargin));
if nargin < 9
    Ns = 200;
end
if nargin < 10
    IsDebug = 0;
end
if nargin < 11
    Gt = [];
end
if nargin < 12
    Cx = [];
end
if nargin < 13
    Cy = [];
end
[T, f, h, y, R, Q, Cx, Cy, N] = FilterInputPreprocess(f, h, y, R, Q, Cx, Cy, x0);
particles = mvnrnd(x0', P0, Ns)';
P = reshape(repmat(P0, 1, Ns), N, N, Ns);
if isempty(Gt)
    Gt = @(particles, w) sum(particles * w', 2); %% MMSE estimation
end
W = ones(1, Ns) / Ns;
x = repmat(x0, 1, T + 1);
for t = 2 : T + 1
    Cx_t = Cx{t - 1};
    Cy_t = Cy{t - 1};
    Rt = R{t - 1};
    Qt = Q{t - 1};
    ft = f{t - 1};
    ht = h{t - 1};
    yt = y{t - 1};
    
    if isvector(particles) && IsDebug     
        [pp, xx] = ksdensity(particles, 'weights', W);
         plot(xx, pp, 'b');
    end
    
    for i = 1 : Ns
        if isempty(Cx_t)
            PP = @(xt) ft('p', xt,particles(:, i), Qt);
        else
            PP = @(xt) ft('p', xt, particles(:, i), Qt, Cx_t);
        end
        [particles(:, i) p P(:,:,i)] = ImportanceFunc(ft, ht, particles(:, i), P(:,:,i), yt, Rt, Qt, Cx_t, Cy_t);

%         particles(:, i) = q('DS', 1); %% Prediction     
        if isempty(Cy_t)
            PL = @(yt) ht('p', yt, particles(:, i), Rt);
        else
            PL = @(yt) ht('p', yt,particles(:, i), Rt, Cy_t);
        end

         W(i) = W(i) * PL(yt) * PP(particles(:, i)) / p;
    end
    W = W / sum(W); %% Update on the newest observation
    
    if isvector(particles) && IsDebug  
        hold on;
        [pp, xx] = ksdensity(particles);
        plot(xx, pp, 'r');
        [pp, xx] = ksdensity(particles, 'Weights', W);
        plot(xx, pp, 'm');
        legend('Prior', 'Prediction', 'Filtering');
        hold off;
        pause;
    end
    
    % Resampling
    index = randsample(1 : Ns, Ns, true, W);
    W = ones(1, Ns) / Ns;
    particles = particles(:, index);
 
    % MCMC
    x(:, t) = Gt(particles, W);
end