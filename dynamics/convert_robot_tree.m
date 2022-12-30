function robot = convert_robot_tree(robot_tree)
% 将matlab urdf生成的robot结构体转换成read_dynamic_file的格式
n = length(robot_tree.Bodies);
dof = n - 1;
M = zeros(4,4,dof);
A = zeros(dof,6);
Tb = zeros(4,4,dof);
mass = zeros(dof, 1);
inertia = zeros(3,3,dof);
ME = eye(4);
for i = 1 : n
    if i < n
        mass(i) = robot_tree.Bodies{i}.Mass;
        Tb(:,:,i) = [eye(3), robot_tree.Bodies{i}.CenterOfMass'; 0 0 0 1];
        I = getInertiaMatrix(robot_tree.Bodies{i}.Inertia);
        inertia(:,:,i) = transform_com_inertia_matrix(I, mass(i), Tb(:,:,i));
        M(:,:,i) = robot_tree.Bodies{i}.Joint.JointToParentTransform * Tb(:,:,i);
        Tb(:,:,i) = tform_inv(Tb(:,:,i));
        axis = Tb(1:3,1:3,i) * robot_tree.Bodies{i}.Joint.JointAxis';
        v = -cross(axis, Tb(1:3,4, i));
        if strcmp(robot_tree.Bodies{i}.Joint.Type, 'revolute') == 1
            A(i,:) = [axis;v];
        else
            A(i,:) = [0;0;0;axis];
        end
        if i > 1
            M(:,:,i) = Tb(:,:,i-1) * M(:,:,i);
        end
    else % end-effector frame
        m_e = robot_tree.Bodies{i}.Mass;
        T_e = [eye(3), robot_tree.Bodies{i}.CenterOfMass'; 0 0 0 1];
        I = getInertiaMatrix(robot_tree.Bodies{i}.Inertia);
        I_e = transform_com_inertia_matrix(I, m_e, T_e);
        ME = Tb(:,:,i-1) * robot_tree.Bodies{i}.Joint.JointToParentTransform;
        T = ME * T_e;
        [I, T] = composite_inertia_matrix(inertia(:,:,i-1), mass(i-1), I_e, m_e, T);
        inertia(:,:,i-1) = I;
        mass(i-1) = mass(i-1) + m_e;
        M(:,:,i-1) =  M(:,:,i-1) * T;
        Tinv = tform_inv(T);
        ME = Tinv * ME;
        A(i-1,:) = adjoint_T(Tinv)*A(i-1,:)';
    end
end
robot.dof = dof; %自由度
robot.mass = mass;%质量
robot.inertia = inertia;% 惯性矩阵
robot.A = A;% screw axis
robot.M = M;% relation at zero position
robot.ME = ME;% end-effector frame
robot.jtMechanics = zeros(dof, 3); % joint damping
robot.gravity = [0, 0, -9.8]; % gravity acceleration in base frame

    function I = getInertiaMatrix(Inertia)
        I = [Inertia(1), Inertia(6), Inertia(5);
            Inertia(6), Inertia(2), Inertia(4);
            Inertia(5), Inertia(4), Inertia(3)];
    end
end