function pinvJ = pinv_J(J, M)
tem = M \ J';
pinvJ = tem * inv(J * tem);