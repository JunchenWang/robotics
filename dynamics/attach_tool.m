function robot = attach_tool(robot, TCP, centerofmass, load, R, I)
% TCP, centerofmass, R are in robot.ME 
mass = robot.mass;
inertia = robot.inertia;
A = robot.A;
M = robot.M;
ME = robot.ME;
[I, T] = composite_inertia_matrix(inertia(:,:,end), mass(end), I, load, ME * [R, centerofmass; 0 0 0 1]);
inertia(:,:,end) = I;
M(:,:,end) =  M(:,:,end) * T;
Tinv = tform_inv(T);
ME = Tinv * ME * TCP;
A(end,:) = adjoint_T(Tinv)*A(end,:)';
robot.A = A;
robot.M = M;
robot.ME = ME;
robot.inertia = inertia;
