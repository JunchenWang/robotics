function frame_err = HandEyeError(A, B, fX)
n = size(A, 2);
fA = zeros(6, n -1);
fB = fA;
for i = 1 : n - 1
    fA(:, i) = ConcatenateFrame(A(:, i + 1), InverseFrame(A(:, i)));
    fB(:, i) = ConcatenateFrame(InverseFrame(B(:, i + 1)), B(:, i));
end

AX = zeros(6, n - 1);
XB = AX;
for i = 1 : n - 1
    AX(:, i) = ConcatenateFrame(fA(:, i), fX);
    XB(:, i) = ConcatenateFrame(ConcatenateFrame(fX, fB(:, i)), InverseFrame(fX));
end
frame_err = FrameDistance(fA, XB);
