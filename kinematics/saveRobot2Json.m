function saveRobot2Json(robot, filename)
% json format for matrix is row major!
n = robot.dof;
inertia = zeros(n, 6);
M = zeros(n, 6);
for i = 1 : n
    inertia(i,:) = [robot.inertia(1,1,i), robot.inertia(2,2,i), robot.inertia(3,3,i),robot.inertia(2,3,i),robot.inertia(1,3,i), robot.inertia(1,2,i)];
    M(i,:) = [logR(robot.M(1:3,1:3,i)), robot.M(1:3,4,i)'];
end
robot.M = M;
robot.inertia = inertia;
fid = fopen(filename,'w'); 
encodedJSON = jsonencode(robot, 'PrettyPrint',true); 
fprintf(fid, encodedJSON); 
fclose(fid);