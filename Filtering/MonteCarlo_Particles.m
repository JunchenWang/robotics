function [particles w] = MonteCarlo_Particles(n, g, h)
if nargin == 2
    h = [];
end
if isempty(h)
    particles = randn(1, n);
else
    particles = h('DS', n);
end
w = zeros(1, n);
for i = 1 : n
    p = particles(:, i);
    w(i) = g(p) / h('v', p);
end
w = w / sum(w);