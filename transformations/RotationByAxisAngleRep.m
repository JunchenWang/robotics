function R = RotationByAxisAngleRep(t)
t = t(1 : 3);
% t = t(:);
t_norm = norm(t);
if 0 == t_norm
    R = eye(3);
    return;
end
R = eye(3) + sin(t_norm) / t_norm * SkewMatrix(t) + (1 - cos(t_norm)) / t_norm.^2 * SkewMatrix(t)^2;