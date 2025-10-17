fW2P = Frame_T_Converter(W2P);
n = length(index);
fBB = zeros(6, n);
N = size(X, 2);
for i = 1 : n
    fP2A = InverseFrame(fA(:, i));
    fBB(:, i) = ConcatenateFrame(fW2P, ConcatenateFrame(fP2A, fX));
    T_ext = Frame_T_Converter(fA(:, i));
    MM = K * T_ext(1:3, :);
    xx = MM * [X(:,:,i);ones(1, N)];
    x(:,:,i) = xx(1:2,:) ./ repmat(xx(3,:), 2, 1);
end

for i = 1 : n
    x(:,:,i) = x(:,:,i) + mvnrnd(zeros(1,2), diag([1 1]), N)' ;
    [kB thetaB] = ExtractKTheta(fBB(:, i));
    tB = fBB(4 : 6, i);
    
    thetaB = thetaB + 0.3 / 180 * pi * randn;
    tB= tB +  0.3 * randn(3, 1);
%     kB = kB + C * (-0.5 + rand(3, 1));
    
    fBB(:, i) = [thetaB * kB; tB];
end

[K_hat fX_hat fW2P_hat err] = SimultaneousCamCalibration(X, x, fBB);
[K_hat2 err fAA] = CamCalibration(X, x);
fX_hat2 = HandEyeCalibration(fAA, fB); 
K - K_hat
K - K_hat2
fX - fX_hat
fX - fX_hat2