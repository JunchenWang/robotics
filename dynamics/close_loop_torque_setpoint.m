function tao = close_loop_torque_setpoint(y, qd, qdd, qe, sum_qe, robot, kp, ki, kd)
n = robot.dof;
Kp = kp * eye(n);
Kd = kd * eye(n);
Ki = ki * eye(n);
b =  qdd +  Kp * qe + Ki * sum_qe +  Kd * (qd - y(n+1:end));
tao = mass_matrix(robot, y(1:n)) * b + velocity_torque(robot, y(1:n), y(n+1:end));