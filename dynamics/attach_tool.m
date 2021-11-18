function [mass, inertia, A, M, ME] = attach_tool(mass, inertia, A, M, ME, TCP, centerofmass, load, R, I)
% TCP, centerofmass, R are in ME 
[I, T] = composite_inertia_matrix(inertia(:,:,end), mass(end), I, load, ME * [R, centerofmass; 0 0 0 1]);
inertia(:,:,end) = I;
M(:,:,end) =  M(:,:,end) * T;
Tinv = tform_inv(T);
ME = Tinv * ME * TCP;
A(end,:) = adjoint_T(Tinv)*A(end,:)';
