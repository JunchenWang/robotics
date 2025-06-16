function test_damping_inverse

robot = readRobotJson('panda_correct.json');
q = rand(7,1);
dq = rand(7,1);
dt = 1e-6;
q1 = q + dq * dt;

p1_F = [0, 0, 0]';
p2_F = [0, 0, 0.5]';

T0 = forward_kin_general(robot, q);
R0 = T0(1:3,1:3);
t0 = T0(1:3,4);
P1_0 = R0 * p1_F + t0;
P2_0 = R0 * p2_F + t0;
rcm = P1_0 + (P2_0 - P1_0) * 0.5;
T1 = forward_kin_general(robot, q1);
R1 = T1(1:3,1:3);
t1 = T1(1:3,4);
P1_1 = R1 * p1_F + t1;
P2_1 = R1 * p2_F + t1;

% x = P1_1 + dot(rcm - P1_1, P2_1 - P1_1) / dot(P2_1 - P1_1, P2_1 - P1_1) * (P2_1 - P1_1); 

[J, ~] = rcm_jacobian(robot, q, dq, p1_F, p2_F, rcm + [0,0,1e-3]');
JJt = J * J';
cond(JJt)
x = 1e-1*rand(3,1);
lambda = 0.1;
q1 = J' * (JJt \ x);
q2 = J' * ((JJt' * JJt + lambda * eye(3)) \ (JJt' * x));
q3 = J' * ((JJt + lambda * eye(3)) \ x);
disp(norm(q1));
disp(norm(q2));
disp(norm(q3));



