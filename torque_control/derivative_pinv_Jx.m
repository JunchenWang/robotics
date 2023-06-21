function dpinvJx = derivative_pinv_Jx(J, M, dJ, dM, x)
tem = M \ J';
dtem = -(M \ dM) * (M \ J') + (M \ dJ');
tem2 = J * tem;
% tem2 = inv(J * tem);
% dtem2 = -((J * tem) \ (dJ * tem + J * dtem)) * tem2;
dpinvJx = (dtem  - tem * (tem2 \ (dJ * tem + J * dtem))) * (tem2 \ x);