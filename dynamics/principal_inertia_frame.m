function [Rbc, Ic] = principal_inertia_frame(Ib)
% Ic = Rbc'*Ib*Rbc
% Ic is diagonalized
[Rbc,Ic] = eig(Ib);
[d,ind] = sort(diag(Ic), 'descend');
Rbc = Rbc(:, ind);
Ic = diag(d);