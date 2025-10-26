function MonteCarlo_Test
g = @(x) normpdf(x, 0.5);
n = 1000;
[particles w] = MonteCarlo_Particles(n, g, @h);
% index = randsample(1 : n, n, true, w);
% particles = particles(:, index);
plot(particles, w, 'r*');
x = -3 : 0.01 : 3;
hold on;
plot(x, g(x) / 170, '*');
display(sum(w .* particles));
function p = h(op, x)
if strcmp(op, 'DS')
    p = -3 + 6 * rand(1, x);
else
    p = 1 / 6;
end

% w = zeros(1, n);
% x = randn(1, n);
% for i = 1 : n
%     w(i) = normpdf(x(i));
% end
% x = randsample(x, n, true, w);
% hist(x, 20);
