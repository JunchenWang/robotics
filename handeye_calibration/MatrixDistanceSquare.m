function [dr dt] = MatrixDistanceSquare(f1, f2)
f = f1 - f2;%ConcatenateFrame(InverseFrame(f1), f2);
dr = sum(f(1:3).^2);
dt = sum(f(4:6).^2);