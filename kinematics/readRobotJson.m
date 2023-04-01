function robot = readRobotJson(filename)
fid = fopen(filename,'r'); 
str = '';
line_in = fgets(fid);
while line_in ~= -1
    str = strcat(str, line_in);
    line_in = fgets(fid);
end
fclose(fid);
robot = jsondecode(str);
n = robot.dof;
inertia = zeros(3,3, n);
M = zeros(4,4, n);
for i = 1 : n
    M(:,:,i) = [exp_w(robot.M(i,1:3)), robot.M(i,4:end)'; 0 0 0 1];
    inertia(:,:,i) = [robot.inertia(i, 1), robot.inertia(i,6), robot.inertia(i,5);
            robot.inertia(i,6), robot.inertia(i,2), robot.inertia(i,4);
            robot.inertia(i,5), robot.inertia(i,4), robot.inertia(i,3)];
end
robot.M = M;
robot.inertia = inertia;

