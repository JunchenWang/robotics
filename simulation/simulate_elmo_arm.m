function simulate_elmo_arm
% 双环控制7轴机器人轨迹，状态空间4个变量：位置，速度，速度误差积分，位置误差积分
port = udpport("byte");
% q = query_joints(port);
robot = convert_robot_tree2(importrobot('urdf\elmourdf13\elmourdf13\urdf\elmourdf13.urdf'));
n = robot.dof;
ptp_move(port,[ -0.2799    0.5257    1.1610   -0.2749    2.1135   -0.6466]);
T_Ftcp=eye(4);
T_Ftcp(3,4)=0.15;
robot.TCP = T_Ftcp;

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