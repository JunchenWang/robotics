function robot = convert_robot_tree(robot_tree, flange)
% 将matlab urdf生成的robot结构体转换成read_dynamic_file的格式
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
Tb = zeros(4,4,n);
mass = zeros(n, 1);
inertia = zeros(3,3,n);
com = zeros(n, 3);
ME = eye(4);
Base = eye(4);
for j = 1 : n
    i = bodylist{j};
    JointTransform = Base * robot_tree.Bodies{i}.Joint.JointToParentTransform;
    if strcmp(robot_tree.Bodies{i}.Joint.Type, 'fixed') == 1
        Base = Base * robot_tree.Bodies{i}.Joint.JointToParentTransform;
        if dof > 0
            m_e = robot_tree.Bodies{i}.Mass;
            T_e = [eye(3), robot_tree.Bodies{i}.CenterOfMass'; 0 0 0 1];
            I = getInertiaMatrix(robot_tree.Bodies{i}.Inertia);
            % matlab读入urdf时，将惯量矩阵转到了link坐标系下，而非质心坐标系，所以要转换一下
            % urdf里的惯量矩阵是关于质心坐标系的
            I_e = transform_com_inertia_matrix(I, m_e, T_e);
            ME = Tb(:,:,dof) * JointTransform;
            T = ME * T_e;
            [I, T] = composite_inertia_matrix(inertia(:,:,dof), mass(dof), I_e, m_e, T);
            inertia(:,:,dof) = I;
            mass(dof) = mass(dof) + m_e;
            com(dof,:) = com(dof,:) + T(1:3,4)';
            M(:,:,dof) =  M(:,:,dof) * T;
            Tinv = tform_inv(T);
            Tb(:,:,dof) = Tinv * Tb(:,:,dof);
            ME = Tinv * ME;
            A(dof,:) = adjoint_T(Tinv)*A(dof,:)';
        end
    else
        dof = dof + 1;
        mass(dof) = robot_tree.Bodies{i}.Mass;
        com(dof,:) = robot_tree.Bodies{i}.CenterOfMass;
        Tb(:,:,dof) = [eye(3), robot_tree.Bodies{i}.CenterOfMass'; 0 0 0 1];
        I = getInertiaMatrix(robot_tree.Bodies{i}.Inertia);
        % matlab读入urdf时，将惯量矩阵转到了link坐标系下，而非质心坐标系，所以要转换一下
        % urdf里的惯量矩阵是关于质心坐标系的
        inertia(:,:,dof) = transform_com_inertia_matrix(I, mass(dof), Tb(:,:,dof));
        M(:,:,dof) = JointTransform * Tb(:,:,dof);
        Tb(:,:,dof) = tform_inv(Tb(:,:,dof));
        axis = Tb(1:3,1:3,dof) * robot_tree.Bodies{i}.Joint.JointAxis';
        v = -cross(axis, Tb(1:3,4, dof));
        if strcmp(robot_tree.Bodies{i}.Joint.Type, 'revolute') == 1
            A(dof,:) = [axis;v];
        else
            A(dof,:) = [0;0;0;axis];
        end
        if dof > 1
            M(:,:,dof) = Tb(:,:,dof-1) * M(:,:,dof);
        end
        Base = eye(4);
        ME = Tb(:,:,dof);
    end
end
robot.dof = dof; %自由度
robot.mass = mass(1:dof);%质量
robot.inertia = inertia(:,:,1:dof);% 惯性矩阵
robot.A = A(1:dof,:);% screw axis
robot.M = M(:,:,1:dof);% relation at zero position
robot.ME = ME;% end-effector frame
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
