function fX = AXXB2(fA, fB)
n = size(fA, 2);
fX_samples = zeros(6, n - 1);
for i = 2 : n
    fA1 = fA(:, i - 1);
    fA2 = fA(:, i);
    fB1 = fB(:, i - 1);
    fB2 = fB(:, i);
    [kA1 thetaA1] = ExtractKTheta(fA1); %#ok<NASGU>
    [kB1 thetaB1] = ExtractKTheta(fB1); %#ok<NASGU>
    [kA2 thetaA2] = ExtractKTheta(fA2); %#ok<NASGU>
    [kB2 thetaB2] = ExtractKTheta(fB2); %#ok<NASGU>
    if norm(cross(kA1, kA2)) < 1e-3
        warning('kA1 and kA2 are approximately parallel');
    end
    kB3 = cross(kB1, kB2);
    kB3 = kB3 / norm(kB3);
    kB2 = cross(kB3, kB1);
    kA3 = cross(kA1, kA2);
    kA3 = kA3 / norm(kA3);
    kA2 = cross(kA3, kA1);
    rX = AngleAxisFromRotation([kA1 kA2 kA3] * [kB1 kB2 kB3]');
    tX = Estimate_tX(fA, fB, rX);
    fX_samples(:, i - 1) = [rX; tX];
end
fX = mean(fX_samples, 2);