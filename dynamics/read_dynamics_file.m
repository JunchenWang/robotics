function robot = read_dynamics_file(filename)
fid = fopen(filename);
n = fscanf(fid, '%f', 1);
mass = zeros(n, 1);
inertia = zeros(3,3,n);
A = zeros(n, 6);
M = zeros(4,4,n);
damp = zeros(n, 1);
for i = 1 : n
    inert = fscanf(fid, '%f %f %f', [7, 1])';
    mass(i) = inert(1);
    inertia(:,:,i) = [inert(2), inert(3), inert(4);
        inert(3), inert(5), inert(6);
        inert(4), inert(6), inert(7)];
    qt = fscanf(fid, '%f %f %f %f %f %f', [7, 1])';
    M(:,:,i) = make_tform_qt(qt);
    A(i, :) = fscanf(fid, '%f %f %f %f %f %f', [6, 1])';
    damp(i) = fscanf(fid, '%f', 1);
end
inert = fscanf(fid, '%f %f %f', [7, 1])';
if ~isempty(inert)
    m_e = inert(1);
    I_e = [inert(2), inert(3), inert(4);
        inert(3), inert(5), inert(6);
        inert(4), inert(6), inert(7)];
    T = make_tform_qt(fscanf(fid, '%f %f %f %f %f %f', [7, 1])');
    ME = make_tform_qt(fscanf(fid, '%f %f %f %f %f %f', [7, 1])');
    [I, T] = composite_inertia_matrix(inertia(:,:,n), mass(n), I_e, m_e, T);
    inertia(:,:,n) = I;
    M(:,:,n) =  M(:,:,n) * T;
    Tinv = tform_inv(T);
    ME = Tinv * ME;
    A(n,:) = adjoint_T(Tinv)*A(n,:)';
else
    ME = eye(4);
end
fclose(fid);
robot.dof = n; %自由度
robot.mass = mass;%质量
robot.inertia = inertia;% 惯性矩阵
robot.A = A;% screw axis
robot.M = M;% relation at zero position
robot.ME = ME;% end-effector frame
robot.damp = damp; % joint damping
robot.gravity = [0, 0, -9.8]; % gravity acceleration in base frame
