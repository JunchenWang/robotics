function [fX, frame_err] = HandEyeCalibration(A, B, flag)
if nargin < 3
    flag = 1;
end
n = size(A, 2);
fA = zeros(6, n -1);
fB = fA;
for i = 1 : n - 1
%     fA(:, i) = ConcatenateFrame(InverseFrame(A(:, i)), A(:, i + 1));
    fA(:, i) = ConcatenateFrame(InverseFrame(A(:, i + 1)), A(:, i));
    fB(:, i) = ConcatenateFrame(InverseFrame(B(:, i)), B(:, i + 1));
end
if flag == 1
    fX = NonlinearAXXB(fA, fB);
elseif flag == 2
    fX = AXXB_de(fA, fB);
elseif flag == 3
    fX = AXXB_dq(fA, fB);
elseif flag == 4
    fX = AXXB_UKF(fA, fB);
elseif flag == 5
    fX = AXXB_optimal(fA, fB);
else
    error('bad choice');
end

AX = zeros(6, n - 1);
XB = AX;
for i = 1 : n - 1
    AX(:, i) = ConcatenateFrame(fA(:, i), fX);
    XB(:, i) = ConcatenateFrame(fX, fB(:, i));
end
frame_err = FrameDistance(AX, XB);
% display(frame_err);
