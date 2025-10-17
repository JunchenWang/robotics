function [J, dJ] = j_matrix(robot, q, qd)
L = robot.L;
c1 = cos(q(1));
c12 = cos(q(1) + q(2));
s1 = sin(q(1));
s12 = sin(q(1) + q(2));
J = [-L * (c1 + c12), -L * c12; -L * (s1 + s12), -L * s12];

dJ = [L * (s1 * qd(1) + s12 * (qd(1) + qd(2))), L * s12 * (qd(1) + qd(2));...
      -L * (c1 * qd(1) + c12 * (qd(1) + qd(2))), -L * c12 * (qd(1) + qd(2))];