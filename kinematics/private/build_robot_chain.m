function [robot, joint_names] = build_robot_chain(body_list)
%BUILD_ROBOT_CHAIN Convert an ordered list of bodies into the robot struct.
n = length(body_list);
dof = 0;
M = zeros(4,4,n);
A = zeros(n,6);
mass = zeros(n, 1);
inertia = zeros(3,3,n);
com = zeros(n, 3);
Base = eye(4);
joint_names = cell(1, n);
for i = 1 : n
    link = body_list{i};
    JointTransform = Base * link.Joint.JointToParentTransform * link.Joint.ChildToJointTransform;
    if strcmp(link.Joint.Type, 'fixed') == 1
        Base = JointTransform;
        m2 = link.Mass;
        if dof > 0 && m2 > 0
            Tcom2 = [eye(3), link.CenterOfMass'; 0 0 0 1];
            I2 = getInertiaMatrix(link.Inertia);
            I2 = transform_com_inertia_matrix(I2, m2, Tcom2);
            T = JointTransform * Tcom2;
            inertia(:,:,dof) = inertia(:,:,dof) + transform_inertia_matrix(I2, m2, T);
            rho = m2 / (mass(dof) + m2);
            com(dof,:) = rho * T(1:3,4)' + (1 - rho) * com(dof,:);
            mass(dof) = mass(dof) + m2;
        end
    else
        dof = dof + 1;
        mass(dof) = link.Mass;
        com(dof,:) = link.CenterOfMass;
        inertia(:,:,dof) = getInertiaMatrix(link.Inertia);
        M(:,:,dof) = JointTransform;
        axis = link.Joint.JointAxis';
        if strcmp(link.Joint.Type, 'revolute')
            A(dof,:) = adjoint_T(tform_inv(link.Joint.ChildToJointTransform)) * [axis;0;0;0];
        elseif strcmp(link.Joint.Type, 'prismatic')
            A(dof,:) = adjoint_T(tform_inv(link.Joint.ChildToJointTransform)) * [0;0;0;axis];
        else
            error("wrong joint type");
        end
        Base = eye(4);
        joint_names{dof} = link.Joint.Name;
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
joint_names = joint_names(1:dof);
end

%% getInertiaMatrix
function I = getInertiaMatrix(Inertia)
        I = [Inertia(1), Inertia(6), Inertia(5);
            Inertia(6), Inertia(2), Inertia(4);
            Inertia(5), Inertia(4), Inertia(3)];
end