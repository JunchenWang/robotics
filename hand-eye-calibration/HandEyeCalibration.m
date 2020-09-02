function [fX, frame_err] = HandEyeCalibration(A, B, flag)
if nargin < 3
    flag = 5;
end
[fA, fB] = selection_process(A, B);
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
n = size(fA, 2);
AX = zeros(6, n);
XB = AX;
for i = 1 : n
    AX(:, i) = ConcatenateFrame(fA(:, i), fX);
    XB(:, i) = ConcatenateFrame(fX, fB(:, i));
end
frame_err = FrameDistance(AX, XB);
% display(frame_err);

function [fA, fB] = selection_process(A, B)
A = EulZYXFrame2RVFrame(A);
B = EulZYXFrame2RVFrame(B);
n = size(A, 2);
fA = zeros(6, n * n);
fB = zeros(6, n * n);
cnt = 1;
for i = 1 : n
    for j = i + 1 : n
        fa = ConcatenateFrame(InverseFrame(A(:, i)), A(:, j));
        fb = ConcatenateFrame(InverseFrame(B(:, i)), B(:, j));
        anga = abs(fa(1:3)) / pi * 180;
        angb = abs(fb(1:3)) / pi * 180;
        if abs(anga - angb) < 1 && abs(anga) > 20 && abs(anga) < 100
            fA(:,cnt) = fa;
            fB(:,cnt) = fb;
            cnt = cnt + 1;
        end
    end
end
fA(:,cnt + 1:end) = [];
fB(:,cnt + 1:end) = [];