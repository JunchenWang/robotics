function d = DistanceR(rA, rB)
if ~isvector(rA)
    rA = AngleAxisFromRotation(rA);
end

if ~isvector(rB)
    rB = AngleAxisFromRotation(rB);
end

[kA thetaA] = ExtractKTheta(rA);
[kB thetaB] = ExtractKTheta(rB);
d = [0 0]';
d(1) = atan2(norm(cross(kA, kB)), dot(kA, kB)) / pi * 180;
d(2) = abs(thetaA - thetaB) / pi * 180;
