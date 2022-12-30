function rpy = R2RPY(R)
theta = R2EulZYX(R);
rpy = [theta(3);theta(2); theta(1)];