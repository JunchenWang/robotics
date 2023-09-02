function robot = convert_robot_tree2(robot_tree, flange)
% 动力学与运动学参数相互分离！
n = robot_tree.NumBodies;
if nargin < 2
    flange_no = n;
else
    flange_no = 0;
    for i = 1 : n
        if strcmp(robot_tree.Bodies{i}.Name, flange) == 1
            flange_no = i;
            break;
        end
    end
    if ~flange_no
        error('no flange name found');
    end
end
i = flange_no;
bodylist = {};
while i
    bodylist = [{i}, bodylist];
    i = findBody(robot_tree.Bodies{i}.Parent);
end
n = length(bodylist);
dof = 0;
M = zeros(4,4,n);
A = zeros(n,6);
mass = zeros(n, 1);
inertia = zeros(3,3,n);
com = zeros(n, 3);
Base = eye(4);
for j = 1 : n
    i = bodylist{j};
    JointTransform = Base * robot_tree.Bodies{i}.Joint.JointToParentTransform;
    if strcmp(robot_tree.Bodies{i}.Joint.Type, 'fixed') == 1
        Base = JointTransform;
        m2 = robot_tree.Bodies{i}.Mass;
        if dof > 0 && m2 > 0
            Tcom2 = [eye(3), robot_tree.Bodies{i}.CenterOfMass'; 0 0 0 1];
            I2 = getInertiaMatrix(robot_tree.Bodies{i}.Inertia);
            I2 = transform_com_inertia_matrix(I2, m2, Tcom2);
            T = JointTransform * Tcom2;
            inertia(:,:,dof) = inertia(:,:,dof) + transform_inertia_matrix(I2, m2, T);
            rho = m2 / (mass(dof) + m2);
            com(dof,:) = rho * T(1:3,4)' + (1 - rho) * com(dof,:);
            mass(dof) = mass(dof) + m2;
        end
    else
        dof = dof + 1;
        mass(dof) = robot_tree.Bodies{i}.Mass;
        com(dof,:) = robot_tree.Bodies{i}.CenterOfMass;
        inertia(:,:,dof) = getInertiaMatrix(robot_tree.Bodies{i}.Inertia);
        M(:,:,dof) = JointTransform;
        axis = robot_tree.Bodies{i}.Joint.JointAxis';
        if strcmp(robot_tree.Bodies{i}.Joint.Type, 'revolute') == 1
            A(dof,:) = [axis;0;0;0];
        else
            A(dof,:) = [0;0;0;axis];
        end
        Base = eye(4);
    end
end
robot.dof = dof; %自由度
robot.mass = mass(1:dof);%质量
robot.inertia = inertia(:,:,1:dof);% 惯性矩阵
robot.A = A(1:dof,:);% screw axis
robot.M = M(:,:,1:dof);% relation at zero position
robot.ME = Base;% end-effector frame
robot.com = com(1:dof,:); % center of mass
robot.gravity = [0, 0, -9.8]; % gravity acceleration in base frame
robot.TCP = eye(4);
    function I = getInertiaMatrix(Inertia)
        I = [Inertia(1), Inertia(6), Inertia(5);
            Inertia(6), Inertia(2), Inertia(4);
            Inertia(5), Inertia(4), Inertia(3)];
    end
    function ind = findBody(obj)
        for k = 1 : n
            if robot_tree.Bodies{k} == obj
                ind = k;
                return;
            end
        end
        ind = 0;
    end
end
