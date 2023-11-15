function simulate_task_control_cuttingline2
% 双环控制7轴机器人轨迹，状态空间4个变量：位置，速度，速度误差积分，位置误差积分
port = udpport("byte");
robot = convert_robot_tree2(importrobot('urdf\ur_5e-calibrated\ur_description\urdf\ur5e-A302.urdf'));
n = robot.dof;
% param = [0.425, 0.3922, 0.1625, 0.1333, 0.0997, 0.0996];
ptp_move(port,[-0.0479   -2.2554   -1.3722   -1.0293    1.5965    0.7480]);
T_Ftcp=eye(4);
T_Ftcp(3,4)=0.15;
robot.TCP = T_Ftcp;
opts = odeset('OutputFcn', @(t, y, flag) odeplot(t, y, flag, port, robot));
y0 = zeros(3 * n + 6, 1);
Kp_s = [100,100,100,100,100,100]';
Ki_s = [50,50,50,100,100,100]';

Kp_p = 5 * [1, 1, 1, 1, 1, 1]';
Ki_p = 0 * [1,1,1,1,1,1]';
Kd_p = 0 * [1,1,1,1,1,1]';

d = @(t, y) [6;5;4;3;2;1]*10;
J = [8;5;4;3;2;2];
B = 2;
r = 200;
tspan = [0, 1000];
 %% 对姿态插值（测试）
% 
% path_locs = load("polyline-surfacetest02.txt");
% path_norms = load("polylinenorm-surfacetest02.txt");
% 
% path_interp = load("polyline-surfacetest-spline025.txt");
% path_cnt = size(path_locs, 1);
% path_cnt_interp = size(path_interp,1);
% pose_matrices = cell(1,path_cnt);
% pose_matrices_interp = cell(1,path_cnt_interp);
% 
% for i = 1:path_cnt
%     z_axis = path_norms(i,:);
%     z_axis = z_axis/norm(z_axis);
%     if i~=path_cnt
%         x_reference = path_locs(i+1,:)-path_locs(i,:);
%     else
%         x_reference = path_locs(i,:)-path_locs(i-1,:);
%     end
%         x_reference = x_reference/norm(x_reference);
%         y_axis = cross(z_axis,x_reference);
%         y_axis = y_axis/norm(y_axis);
%         x_axis = cross(y_axis,z_axis);
%         x_axis = x_axis/norm(x_axis);
%         R=[x_axis',y_axis',z_axis'];
%         loc_t=path_locs(i,:).';
%         loc_t=loc_t./1000.0;
%         pose_matrix=[R,loc_t;0,0,0,1];
%         pose_matrices{i} = pose_matrix;    
% end
% 
% % 对姿态插值
% 
% path_cnt_interp = 0;
% R_start0 = pose_matrices{1}(1:3,1:3);
% for i = 1 : path_cnt-1
%     delta_t = 0;
%     R_start = pose_matrices{i}(1:3,1:3);
%     R_end = pose_matrices{i+1}(1:3,1:3);
%     delta_h = path_locs(i+1,1)-path_locs(i,1);
%     while delta_t <= delta_h
%         R_temp = R_start *exp_w(logR(R_start'*R_end)*delta_t/delta_h);
%         delta_t = delta_t+0.2;
%         path_cnt_interp = path_cnt_interp + 1;
%         pose_matrices_interp{path_cnt_interp} = [R_temp,path_interp(path_cnt_interp,:)'./1000.0;0,0,0,1];
%         % pose_matrices_interp{path_cnt_interp} = [R_start0,path_interp(path_cnt_interp,:)'./1000.0;0,0,0,1];
% 
%     end
% end
% pose_matrices = pose_matrices_interp;
%% 直接用C++计算得到的轨迹控制
pose_matrices_file = load("F:\MICR\temp\trajInterp.txt");
path_cnt_interp = size(pose_matrices_file,1)/4;
pose_matrices = cell(1,path_cnt_interp);
for i = 1:path_cnt_interp
    pose_matrices{i}=pose_matrices_file(4*(i-1)+1:4*i,:);
end
T0 = pose_matrices{1};
tol = [1e-5,1e-5];
start = query_joints(port);
[q0, flag1] = inverse_kin_general(robot,T0,start,tol);
% % y0(1:n) = [0 75 0 -94 0 -81 ] / 180 * pi;
ptp_move(port, q0);
% kesai = cal_kuka_kesai(y0);
R = pose_matrices{1}(1:3,1:3);
for i = 2:path_cnt_interp
    Ti = pose_matrices{i};
    Ti(1:3,1:3) = R;
    start = query_joints(port);
    [y0(1:n), flag] = inverse_kin_general(robot,Ti,start,tol);
    if ~flag
        error('no solution');
    end
    Freq = 200;
    rate = rateControl(Freq);
    set_joints(port, y0(1:n));
    waitfor(rate);

end
end

function [u, e] = speed_controller(t, y, desired_speed, Kp, Ki)
% y(1): theta, y(2) : w, y(3) sum e
n = round((length(y) + 1) / 4);
[wd, e2] = desired_speed(t, y);
e1 = wd - y(n + 1: 2 * n);
u = Kp .* e1 + Ki .* y(2 * n + 1 : 3 * n);
e = [e1;e2];
end

function [desired_speed, e] = task_pos_controller(t, y, desired_task_pos, robot, Kp, Ki, Kd)
n = round((length(y) + 1) / 4);
q = y(1:n);
qd = y(n + 1 : 2*n);
Jb = jacobian_matrix(robot, q);
V = Jb * qd;
T = forward_kin_general(robot, q);
R = T(1:3,1:3);
p = T(1:3,4);
[Td, Vd] = desired_task_pos(t);
% Td =Td{1,1}; % 改
wd = Vd(1:3);
vd = Vd(4:6);
Rd = Td(1:3,1:3);
pd = Td(1:3,4);
xe = logR(R'*Rd)';
pe = pd - p;
wb = R' * wd + Kp(1:3) .* xe + Ki(1:3) .* y(3*n + 1 : 3*n + 3) + Kd(1:3) .* (R' * wd - V(1:3));
vb = R'*(vd + Kp(4:6) .* pe + Ki(4:6) .* y(3*n + 4 : 3*n + 6) + Kd(4:6) .* (vd - R * V(4:6)));
e = [xe; pe];
desired_speed = lsqminnorm(Jb,[wb;vb]);
% disp(Jb * pinv(Jb));
end


function [Td, Vd] = desired_task_pos(t, Ts, Te, pp, ppd, ppdd)
s = ppval(pp, t);
sd = ppval(ppd, t);
sdd = ppval(ppdd, t);
ps = Ts(1:3,4);
pe = Te(1:3,4);
Rs = Ts(1:3,1:3);
Re = Te(1:3,1:3);
pd = ps + (pe - ps)*s;
xe = logR(Rs'*Re)';
Rd = Rs * exp_w(xe*s);
Td = [Rd, pd; 0 0 0 1];
Vd = [Rs * xe; (pe - ps)] * sd;
end


% here
function [Td, Vd] = desired_task_pos_cuttingline(t, path)
dt = 0.02;
cnt = round(t / dt) + 1;
Td = path{cnt};
%Td(1:3,1:3) = path{1}(1:3,1:3);
Vd = cartesian_velocity_estimator(t, Td);
end




function ret = odeplot(t, y, flag, port, robot)
if strcmp(flag, 'init') == 1
elseif isempty(flag)
    set_joints(port, y(1:robot.dof, end));
else
end
ret = 0;
end


function yd = joint_motor_dynamic(t, y, u, d, J, B, r)
% y(1): theta, y(2): w, y(3): sum e
n = round((length(y) + 1) / 4);
yd = zeros(3 * n + 6, 1);
[U, e] = u(t,y);
yd(1:n) = y(n +1 : 2 * n);
yd(n +1 : 2*n) = (U - d(t, y) ./ r - B .* y(n +1: 2 *n)) ./ J;
yd(2 * n + 1 : end) = e;
% disp(u(t,y));
end

function dist = dist2line(pt, p1, p2)
    s = p2 - p1;
    s = s / norm(s);
    l1 = pt - p1;
    l2 = dot(l1, s);
    dist = real(sqrt(l1'*l1 - l2*l2));
end