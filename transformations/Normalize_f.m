function fn = Normalize_f(f)
f = reshape(f, 6, []);
N = size(f, 2);
fn = zeros(6, N);
fn(1 : 3, :) = Normalize_r(f(1 : 3, :));
fn(4 : 6, :) = f(4 : 6, :);