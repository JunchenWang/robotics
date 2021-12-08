function [vel, acc] = cartesian_velocity_estimator(t, T)
persistent prev_t;
persistent prev_T;
persistent prev_V;
vel = zeros(6,1);
acc = zeros(6,1);
if isempty(prev_t)
    prev_t = t;
    prev_T = T;
    prev_V = zeros(6,1);
else
    delta = t - prev_t;
    dT = (T - prev_T) / delta;
    se_V = tform_inv(T) * dT;
    V = [skew_mat2vec(se_V(1:3,1:3))';se_V(1:3, 4)];
    dV = (V - prev_V) / delta;
    vel(1:3) = T(1:3,1:3) * V(1:3);
    vel(4:6) = T(1:3,1:3) * V(4:6);
    acc(1:3) = T(1:3,1:3) * dV(1:3);
    acc(4:6) = T(1:3,1:3) * (dV(4:6) + cross(vel(1:3), vel(4:6)));
    prev_t = t;
    prev_T = T;
    prev_V = V;
end

