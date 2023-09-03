function param = get_dynamics(robot)
n = robot.dof;
mass = robot.mass;
com = robot.com;
inertia = robot.inertia;
param = zeros(10 * n, 1);
for i = 1 : n
    param(10 * (i -1) + 1 : 10 * i) = [mass(i), mass(i) * com(i,:), inertia(1,1,i), inertia(2,2,i), inertia(3,3,i),...
                inertia(1,2,i), inertia(1,3,i), inertia(2,3,i)];
end