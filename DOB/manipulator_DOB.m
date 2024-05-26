function [tau, dz] = manipulator_DOB(robot, Y, tau, M, C, G, J, y)
% Y: DOB gain
% y: state
n = robot.dof;
dq = y(n + 1 : 2 * n);
z = y(2 * n + 1 : 3 * n);
P = Y * dq;
tau_d = z + P;
tau_d = J' * ((J * (M \ J'))  \ (J * (M \ tau_d))); % minus the component projected into the null space !
tau = tau - tau_d;
dz = Y * (M \ (C * dq + G - tau - P - z));
