function [err_relative_r, err_relative_t, fX] = testAXXB
n = 200;
[fA, fB, fX] = GenerateABX(n);
fA_n = fA;
fB_n = fB;
MC = 1;
r = norm(fX(1:3));
t = norm(fX(4:6));
err_r = zeros(3, MC);
err_t = err_r;
err_relative_r = zeros(3, 10);
err_relative_t = zeros(3, 10);
for C = 1 : 1
    time = tic;
for k = 1 : MC
    for i = 1 : n
        
        rB = fB(1:3, i);
        tB = fB(4 : 6, i);
        rA = fA(1:3, i);
        tA = fA(4 : 6, i);
        %add noise
        tA= tA +  0.5 * C * 0.2 * randn(3, 1);
        rA = rA +0.5 * C * 0.01 * randn(3, 1);
        tB= tB + 0.5 * C * 0.2 * randn(3, 1);
        rB = rB + 0.5 * C * 0.01 * randn(3, 1);
              
        fA_n(:, i) = [rA; tA];
        fB_n(:, i) = [rB; tB];
    end
    disp('***************optimal*******************');
    fX1 = AXXB_optimal(fA_n, fB_n) - fX;
    disp('***************de*******************');
    fX2 = AXXB_de(fA_n, fB_n) - fX;
    disp('***************nonlinear*******************');
    fX3 = NonlinearAXXB(fA_n, fB_n) - fX;
    disp('***************finish*******************');
    err_r(1,k) = norm(fX1(1:3));
    err_t(1,k) = norm(fX1(4:6));
    err_r(2,k) = norm(fX2(1:3));
    err_t(2,k) = norm(fX2(4:6));
    err_r(3,k) = norm(fX3(1:3));
    err_t(3,k) = norm(fX3(4:6));
end

err_rotation = sqrt(mean(err_r.^2, 2));
err_translation = sqrt(mean(err_t.^2, 2));
err = [err_rotation, err_translation];
display(err);
err_relative_r(:, C) = err_rotation;
err_relative_t(:, C) = err_translation;
toc(time);
end
plot(err_relative_r');
figure;
plot(err_relative_t');
disp(fX);
% fX1 = NonlinearAXXB(fA, fB);
% fX2 = AXXB_de(fA, fB);
% fX3 = AXXB_dq(fA, fB);
% fX4 = AXXB_dq2(fA, fB);
% fX5 = AXXB_optimal(fA, fB);
% fX1
% fX2
% display(fX);
% err_wang = (fX5 - fX)';
% err_no = (fX1 - fX)';
% err_de = (fX2 - fX)';
% err_dq = (fX3 - fX)';
% err_dq2 = (fX4 - fX)';
% display(err_no);
% display(err_de);
% display(err_dq);
% display(err_dq2);
% display(err_wang);



function [fA, fB, fX] = GenerateABX(n)
fA = zeros(6, n);
fB = fA;
% tX = -100 * ones(3, 1) + rand(3, 1) * 200;
tX = [40, 35, -360]';
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



function [err_UKF, err_ls] = test_performance(MC)
 n = 30;
 fid1 = fopen('eyehand_ukf.txt','w');
 fid2 = fopen('eyehand_ls.txt','w');
 
 [fA, fB, fX] = GenerateABX(n);
 sigma_tA = diag([9 9 9]);
 sigma_rA = diag([9 9 9]);
 sigma_tB = diag([1 1 1]);
 sigma_rB = diag([1 1 1]);
 %add noise
 err_UKF = zeros(6, MC);
 err_ls = zeros(6, MC);
 x0 = ones(6, 1);
 P0 = eye(6);
 q = 0.1;
 Q = q * blkdiag(sigma_rB + sigma_rA, sigma_tB + sigma_tA);
 R = cell(1, n);
 R{1} = blkdiag(sigma_rB, sigma_tB, sigma_rA, sigma_tA);
 for i = 2 : n
     R{i} = blkdiag(R{i - 1}, sigma_rB, sigma_tB, sigma_rA, sigma_tA);
 end
 fprintf(fid1, '%f\n %f %f %f %f %f %f\n\n', [q; fX]);
 fprintf(fid2, '%f\n %f %f %f %f %f %f\n\n', [q; fX]);
 try
     for c = 1 : MC
         t = tic;
         fA_n = mvnrnd(fA', blkdiag(sigma_rA, sigma_tA))';
         fB_n = mvnrnd(fB', blkdiag(sigma_rB, sigma_tB))';
         Cy = cell(1, n);
         Cy{1} = fB_n(:, 1);
         y = cell(1, n);
         y{1} = fA_n(:, 1);
         for i = 2 : n
             Cy{i} = [Cy{i - 1}; fB_n(:, i)];
             y{i} = [y{i - 1}; fA_n(:, i)];
         end
         
         x = UKalmanFilter(@f, @h, x0, P0, y, R, Q, [], Cy, [0.25 3 0]);        
         plot(x');
         pause(0.1);
         UKF_x = x(:, end);
         UKF_x = ThetaFrame_T_Converter(ThetaFrame_T_Converter(UKF_x));
         err_UKF(:, c) = UKF_x - fX;
         
         fA_n = ThetaFrame2RVFrame(fA_n);
         fB_n = ThetaFrame2RVFrame(fB_n);
         fX_ls = RVFrame2ThetaFrame(AXXB_LS(fA_n, fB_n));
         err_ls(:, c) = fX_ls - fX;
         fprintf(fid1, '%f %f %f %f %f %f\n', err_UKF(:, c));
         fprintf(fid2, '%f %f %f %f %f %f\n', err_ls(:, c));
         elapsedtime = toc(t);
         estimatedtime = (MC - c) * elapsedtime;
         display(estimatedtime / 60);
     end
 catch me
     display(me.message);
     fclose(fid1);
     fclose(fid2);
 end





