function [robot, joint_names] = get_robot_chain(robot_tree, flange_link, base_link)
% 创建从base_link 到 flange_link的运动链
robot = [];
joint_names = {};
%% find base and flange
if nargin < 3 || strcmp(base_link, robot_tree.BaseName)
    base = robot_tree.Base;
else
    base_id = find(strcmp(robot_tree.BodyNames, base_link));
    if isempty(base_id)
        disp('not a valid base');
        return;
    else
        base = robot_tree.Bodies{base_id};
    end
end

if nargin < 2 
    flange = robot_tree.Bodies{end};
else
    flange_id = find(strcmp(robot_tree.BodyNames, flange_link));
    if isempty(flange_id)
        disp('not a valid flange');
        return;
    else
        flange = robot_tree.Bodies{flange_id};
    end
end

%% check chain
link = flange;
flag = 0;
maxN = length(robot_tree.Bodies);
body_list = cell(1, maxN);
head = maxN + 1;
while(~isempty(link))
    if(link == base)
        flag = 1;
        break;
    end 
    head = head - 1;
    body_list{head} = link;
    link = link.Parent;
end
body_list = body_list(head:end);
n = length(body_list);
if ~flag || n == 0
    disp('not a valid chain');
    return;
end

%% get parameters
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