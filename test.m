function test(N)
if nargin < 1
    N = 10;
end
%% twist test
disp('twist test');
v = zeros(N, 6);
for i = 1 : N - 1
w = rand(1,3);
v(i, :) = [w / norm(w) * rand * pi, rand(1, 3) * 1000];
end
T = exp_twist(v);
vv = logT(T);
disp(norm(vv - v, 'fro') / sqrt(N));

%% jacobian test
disp('jacobian test');
L = 0.5;
B1 = [0 0 1 0 2*L 0];
B2 = [0 0 1 0 L 0];
B3 = [0 0 1 0 0 0];
B = [B1; B2; B3];
S1 = [0 0 1 0 0 0];
S2 = [0 0 1 0 -L 0];
S3 = [0 0 1 0 -2*L 0];
S = [S1; S2; S3];
M = make_tform(eye(3), [2 * L, 0, 0]);
config = rand(1, 3);
config = pi * config / norm(config);
T = forward_kin(M, B, config, 'b');
TT = forward_kin(M, S, config, 's');
Jb = jacobian(B, config, 'b');
Js = jacobian(S, config, 's');
disp(norm(Js - adjoint_T(T) * Jb, 'fro'));
disp(norm(TT - T, 'fro'));

%% singularity
config = [0 0 0];
T = forward_kin(M, B, config, 'b');
J = jacobian(B, config, 'b');
J = J([4:5,3],:);
disp(null(J'));