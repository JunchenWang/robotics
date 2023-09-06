function dpinvJ = derivative_pinv_J(J, M, dJ, dM)
tem = M \ J';
dtem = -(M \ dM) * (M \ J') + (M \ dJ');
tem2 = inv(J * tem);
% dtem2 = -((J * tem) \ (dJ * tem + J * dtem)) * tem2;
dpinvJ = (dtem  - tem * ((J * tem) \ (dJ * tem + J * dtem))) * tem2;