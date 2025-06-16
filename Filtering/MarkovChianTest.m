function MarkovChianTest
N = 1e3;
i = 2;
x = zeros(1, N + 1);
while i <= N + 1
  x(i) = x(i - 1) + randn;  
  i = i + 1;
end
plot(x);