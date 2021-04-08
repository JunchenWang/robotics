function kesai = adjust_kesai(a, b, kesai)
K = 0.5;
alpha = 25;
kesai = kesai + K * (b - a) / 2 * (exp(-alpha*(kesai-a)/(b-a))-exp(-alpha*(b-kesai)/(b-a)));