function pinvZ = pinv_Z(Z, M)
pinvZ = (Z' * M * Z) \ (Z' * M);