function fX = AXXB_de(fA, fB)
n = size(fA, 2);
kA = zeros(3, n);
kB = kA;
for i = 1 : n
    kA(:, i) = ExtractKTheta(fA(:, i));
    kB(:, i) = ExtractKTheta(fB(:, i));
end
R = RotationProblem(kA, kB);
rX = AngleAxisFromRotation(R);
tX = Estimate_tX(fA, fB, rX);
fX = [rX; tX];


function tX = Estimate_tX(fA, fB, rX)
n = size(fA, 2);
A = zeros(3 * n, 3);
b = zeros(3 * n, 1);
RX = RotationByAxisAngleRep(rX);
for i = 1 : n
    A(3 * i - 2 : 3 * i, :) = eye(3) - RotationByAxisAngleRep(fA(1 : 3, i));
    b(3 * i - 2 : 3 * i, :) = fA(4 : 6, i) - RX * (fB(4 : 6, i));
end
tX = A \ b;




