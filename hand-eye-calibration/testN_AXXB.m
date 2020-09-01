function [err_relative_r err_relative_t fX] = testN_AXXB(n1, n2)
n = n2;
[fA fB fX] = GenerateABX(n);
MC = 500;
r = norm(fX(1:3));
t = norm(fX(4:6));
err_r = zeros(3, MC);
err_t = err_r;
err_relative_r = zeros(3, n2 - n1 + 1);
err_relative_t = zeros(3, n2 - n1 + 1);
C = 5;
for N = n1 :  n2 
    fA_n = zeros(6, N);
    fB_n = zeros(6, N);
    for k = 1 : MC
        time = tic;
        for i = 1 : N
            rB = fB(1:3, i);
            tB = fB(4 : 6, i);
            rA = fA(1:3, i);
            tA = fA(4 : 6, i);
            %add noise
            tA= tA +  C * 0.2 * randn(3, 1);
            rA = rA + C * 0.01 * randn(3, 1);
            tB= tB + C * 0.2 * randn(3, 1);
            rB = rB + C * 0.01 * randn(3, 1);
            
            fA_n(:, i) = [rA; tA];
            fB_n(:, i) = [rB; tB];
        end
        fX1 = AXXB_optimal(fA_n, fB_n) - fX;
        fX2 = AXXB_de(fA_n, fB_n) - fX;
        fX3 = NonlinearAXXB(fA_n, fB_n) - fX;
        err_r(1,k) = norm(fX1(1:3)) / r;
        err_t(1,k) = norm(fX1(4:6)) / t;
        err_r(2,k) = norm(fX2(1:3)) / r;
        err_t(2,k) = norm(fX2(4:6)) / t;
        err_r(3,k) = norm(fX3(1:3)) / r;
        err_t(3,k) = norm(fX3(4:6)) / t;
        time = toc(time);
        time = (MC - k) * time * (n2 - N + 1) / 3600;
        display(time);
    end
    err_rotation = sqrt(mean(err_r.^2, 2));
    err_translation = sqrt(mean(err_t.^2, 2));
    err = [err_rotation, err_translation];
    display(err);
    err_relative_r(:, N-n1+1) = err_rotation;
    err_relative_t(:, N-n1+1) = err_translation;
    subplot(1,2,1); plot(err_relative_r');
    subplot(1,2,2); plot(err_relative_t');
    
end
% toc(time);


function [fA fB fX] = GenerateABX(n)
fA = zeros(6, n);
fB = fA;
% tX = -100 * ones(3, 1) + rand(3, 1) * 200;
tX = [40, 35, -260]';
kX = rand(3, 1);
kX = kX / norm(kX);
thetaX = rand * pi;
rX = thetaX * kX;
fX = [rX; tX];
for i = 1 : n
    % true value
    tA =[-50 -50 200]' + diag([100, 100, 100]) * rand(3, 1);
    kA = rand(3, 1);
    kA = kA / norm(kA);
    thetaA = rand * pi;
    rA = thetaA * kA;
    fA(:, i) = [rA; tA];
    % calculate B without noise
    fB(:, i) = ConcatenateFrame(ConcatenateFrame(InverseFrame(fX), fA(:, i)), fX);

end




