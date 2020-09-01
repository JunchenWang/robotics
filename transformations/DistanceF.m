function D = DistanceF(fA, fB)
d = DistanceR(fA(1 : 3), fB(1 : 3));
trans = norm(fA(4 : 6) - fB(4 : 6));
D = [d; trans];