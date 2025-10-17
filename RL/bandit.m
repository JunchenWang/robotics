function bandit

k = 10;
run = 2000;
epsion = 0;
steps = 2000;
avg_r = zeros(run, steps);
rng("default")
for j = 1 :run
    Q = zeros(k, 1);
    N = zeros(k, 1);
    E = randn(k, 1);
    for i = 1 : steps
        if rand < epsion
            A = randsample(k, 1);
        else
            [~, A] = max(Q);
        end
        N(A) = N(A) + 1;
        R = randn + E(A);
        Q(A) = Q(A) + 1 / N(A) * (R - Q(A));
        avg_r(j, i) = R;
    end
end
plot(sum(avg_r, 1) / run, 'r', 'LineWidth',1);

rng("default")
epsion = 0.01;
for j = 1 :run
    Q = zeros(k, 1);
    N = zeros(k, 1);
    E = randn(k, 1);
    for i = 1 : steps
        if rand < epsion
            A = randsample(k, 1);
        else
            [~, A] = max(Q);
        end
        N(A) = N(A) + 1;
        R = randn + E(A);
        Q(A) = Q(A) + 1 / N(A) * (R - Q(A));
        avg_r(j, i) = R;
    end
end
hold on;
plot(sum(avg_r, 1) / run, 'g', 'LineWidth',1);
hold off;

rng("default")
epsion = 0.1;
for j = 1 :run
    Q = zeros(k, 1);
    N = zeros(k, 1);
    E = randn(k, 1);
    for i = 1 : steps
        if rand < epsion
            A = randsample(k, 1);
        else
            [~, A] = max(Q);
        end
        N(A) = N(A) + 1;
        R = randn + E(A);
        Q(A) = Q(A) + 1 / N(A) * (R - Q(A));
        avg_r(j, i) = R;
    end
end
hold on;
plot(sum(avg_r, 1) / run, 'b', 'LineWidth',1);
hold off;