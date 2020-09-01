function P_new = ApplyRotationToPoint(r, P, type)
if nargin == 2
    type = 'rv';
end
P = P(:);
if strcmp(type, 'rv')
    phi = norm(r);
    if phi == 0
        P_new = P;
        return;
    end
    crosst_P = cross(r, P);
    P_new = P + sin(phi) / phi * crosst_P + (1 - cos(phi)) / phi.^2 * cross(r, crosst_P);
elseif strcmp(type, 'angle')
    R = thetaXYZ2R(r);
    P_new = R * P;
else
    error('unsupported rotation representation');
end