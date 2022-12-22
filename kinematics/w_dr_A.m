function A = w_dr_A(r)
% wb = A(r)dr/dt
norm_r = norm(r);
if norm_r ~= 0
   sr = vec2skew_mat(r);
    A = eye(3) - (1 - cos(norm_r)) / norm_r^2 * sr + (norm_r - sin(norm_r))/norm_r^3 * sr*sr;
else
    A = eye(3);
end

