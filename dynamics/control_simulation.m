function control_simulation
u = udpport("byte");
robot = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
n = robot.dof;
dynamics = @(t, y, tao) manipulator_dynamics(robot, tao, @ext_wrench, t, y);
MassMatrix = @(t, y) [eye(n), zeros(n); zeros(n), mass_matrix(robot, y(1:n))];

opts = odeset('Mass',MassMatrix,'OutputFcn',@odeplot);
y0 = zeros(2*n,1);
% y0(1:7) = [0;75;0;-90;0;-90;0] / 180 * pi;
ptp(y0(1:n)');
clear computed_torque_control_law;
targetJointPosition  = [pi/2 pi/3 pi/6 2*pi/3 -pi/2 -pi/3, 0]';
targetJointVelocity  = zeros(n,1);
targetJointAcceleration = zeros(n,1);
control_target = @(t, y)computed_torque_control_law(dynamics, robot,...
    @(t) targetJointPosition, @(t) targetJointVelocity,@(t) targetJointAcceleration, t, y);

tspan = [0, 10];
% [t,y] = ode45(@(t, y)dynamics(t, y, @joint_torque),tspan,y0,opts);
[t,y] = ode45(control_target,tspan,y0,opts);
figure;
plot(t, y(:,1:n)');
    function Fext = ext_wrench(t, y)
        Fext = zeros(6,1);
    end
    function ret = odeplot(t, y, flag)
        if strcmp(flag, 'init') == 1
        elseif isempty(flag)
            setJoints(y(1:n));       
        else
        end
        ret = 0;
    end
    function setJoints(jt)
        cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
            ,jt(6), jt(7));
        writeline(u,cmd,"192.168.3.34",7755);
    end
    function ptp(jts, vel)
        if nargin < 2
            vel = .1;
        end
        start = queryJoints;
        init_y = [start'; zeros(n,1)];
        
        clear computed_torque_control_law;
        wayPoints = [start',jts'];
        Freq = 500;
        T = max(abs(jts - start) / vel);
        if T == 0
            return;
        end
        numSamples = max(round(T * Freq) + 1, 2);
        [q,qd,qdd,tSamples,~] = trapveltraj(wayPoints,numSamples,'EndTime',T);
%         fid = fopen('q.txt', 'w');
%         fprintf(fid, '%f %f %f %f %f %f %f\n', q);
%         fclose(fid);
%         
%         fid = fopen('qd.txt','w');
%         fprintf(fid, '%f %f %f %f %f %f %f\n', qd);
%         fclose(fid);
%         
%         fid = fopen('qdd.txt','w');
%         fprintf(fid, '%f %f %f %f %f %f %f\n', qdd);
%         fclose(fid);
        target_q = @(t) interp1(tSamples, q', t)';
        target_qd = @(t) interp1(tSamples, qd', t)';
        target_qdd = @(t) interp1(tSamples, qdd', t)';
        target = @(t, y)computed_torque_control_law(dynamics, robot,...
            target_q, target_qd, target_qdd, t, y);
        [tt,yy] = ode45(target,[0, T],init_y,opts);
        figure;
        plot(tt, yy(:,1:n)');
    end
    function joints = queryJoints
        % ";" 表示查询关节角
        writeline(u,"robot;","192.168.3.34",7755);
        s = readline(u);
        joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
        joints = mod(joints + pi, 2*pi) - pi;
    end
end
