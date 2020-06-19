function Tv = se_twist(v)
if isvector(v)
    Tv = [vec2skew_mat(v(1:3)), [v(4), v(5), v(6)]'; 0 0 0 0];
else
    n = size(v, 1);
    Tv = zeros(4, 4, n);
    skm = vec2skew_mat(v(:, 1 : 3));
    for i = 1 : n
        Tv(1 : 3, 1 : 4, i) = [skm(:, :, i), v(i, 4 : 6)'];
    end
end