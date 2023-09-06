function dpinvZ = derivative_pinv_Z(Z, M, dZ, dM)
tem1 = Z' * M * Z;
dtem1 = dZ' * M * Z + Z' * dM * Z + Z' * M * dZ;
tem2 = Z' * M;
dtem2 = dZ' * M + Z' * dM;
dpinvZ = -(tem1 \ dtem1) * (tem1 \ tem2) + (tem1 \ dtem2); 