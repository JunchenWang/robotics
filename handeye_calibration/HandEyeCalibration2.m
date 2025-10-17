function [fX, frame_err] = HandEyeCalibration2(M, B)
n = size(M, 3);
AA = zeros(6, n -1);
BB = AA;
for i = 1 : n - 1
    M1 = M(:,:,i);
    [K R C] = Decomposition(M1);
    M1 = K * R * [eye(3), -C];
    M2 =  M(:,:,i + 1);   
    [K R C] = Decomposition(M2);
    M2 = K * R * [eye(3), -C];
    N1 = M1(1:3,1:3);
    n1 = M1(:, 4);
    N2 = M2(1:3,1:3);
    n2 = M2(:, 4);
    B1 = B(:, :, i);
    B2 = B(:, :, i + 1);
    N = N1 \ N2;
    [U D V] = svd(N);
    N = U * eye(3) * V';
    AA(:, i) = Frame_T_Converter([N, N1 \ (n2 - n1); zeros(1, 3), 1]);
    BB(:, i) = Frame_T_Converter(B1 * InvertT(B2));
end

fX = NonlinearAXXB(AA, BB);
AX = zeros(6, n - 1);
XB = AX;
for i = 1 : n - 1
    AX(:, i) = ConcatenateFrame(AA(:, i), fX);
    XB(:, i) = ConcatenateFrame(fX, BB(:, i));
end
frame_err = FrameDistance(AX, XB);
display(frame_err);
% fX = AXXB(AA, BB);
% AX = zeros(6, n);
% XB = AX;
% for i = 1 : n - 2
%     AX(:, i) = ConcatenateFrame(AA(:, i), fX);
%     XB(:, i) = ConcatenateFrame(fX, BB(:, i));
% end
% err = ThetaFrameDistance(AX, XB)
