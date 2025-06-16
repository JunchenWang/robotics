function x = DrawSamples
%(pdf, N, a, b)
% samples = linspace(a, b, 200 * N);
% w = pdf(samples);
% w = w / sum(w);
% x = randsample(samples, N, true, w);
p = rand;
x = roots([-0.25 0 0.75 0.5 - p]);
x = x(3);