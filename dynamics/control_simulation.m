function control_simulation
u = udpport("byte");
robot = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
n = robot.dof;
dynamics = @(t, y, tao) manipulator_dynamics(robot, tao, @ext_wrench, t, y);
MassMatrix = @(t, y) [eye(n), zeros(n); zeros(n), mass_matrix(robot, y(1:n))];

opts = odeset('Mass',MassMatrix,'OutputFcn',@odeplot);
y0 = zeros(2*n,1);
% y0(4) = pi / 3;
% targetJointPosition  = [pi/2 pi/3 pi/6 2*pi/3 -pi/2 -pi/3, 0];
ptp(y0(1:n)');

clear computed_torque_control_law;
targetJointPosition  = [pi/2 pi/3 pi/6 2 * pi/3 -pi/2 -pi/3, 0]';
targetJointVelocity  = zeros(n,1);
targetJointAcceleration = zeros(n,1);
control_target = @(t, y)computed_torque_control_law(dynamics, robot,...
    @(t) targetJointPosition, @(t) targetJointVelocity,@(t) targetJointAcceleration, t, y);

% freq = 50;
% time = 10;
% tspan = linspace(0, time, freq * time);
tspan = [0, 10];
[t,y] = ode45(@(t, y)dynamics(t, y, @joint_torque),tspan,y0,opts);
% [t,y] = ode45(control_target,tspan,y0,opts);
figure;
plot(t, y(:,1:n)');

    function tao = joint_torque(t, y)
        %         gt = gravity_torque(mass, inertia, A, M, y(1:n));
        %         gt = gravity_velocity_torque(mass, inertia, A, M, y(1:n), y(n+1:end));
        %         tao = gt + ones(n,1);
        tao = zeros(n,1);
    end
    function Fext = ext_wrench(t, y)
        Fext = zeros(6,1);
    end

    function ret = odeplot(t, y, flag)
        if strcmp(flag, 'init') == 1
            % %             figure;
            %             clf;
            %             axis([t(1), t(2), -2*pi, 2*pi]);
            %             grid on;
            %             hold on;
            %             cnt = cnt +1;
            %             tracer(:, cnt) = y(1:n);
            %             tick(cnt) = t(1);
        elseif isempty(flag)
            setJoints(y(1:n));
            %             step = length(t);
            %             if cnt + step <= MAX
            %                 tracer(:, cnt + 1 : cnt + step) = y(1:n,:);
            %                 tick(cnt +1 : cnt + step) = t;
            %                 cnt = cnt + step;
            %                 plot(tick(1:cnt), tracer(:, 1:cnt));
            %             end
            
        else
            %             hold off;
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
            vel = .4;
        end
        start = queryJoints;
        init_y = [start'; zeros(n,1)];
        
        clear computed_torque_control_law;
        wayPoints = [start',jts'];
        Freq = 200;
        %         r = rateControl(Freq);
        T = max(abs(jts - start) / vel);
        if T == 0
            return;
        end
        numSamples = round(T * Freq) + 2;
        [q,qd,qdd,tSamples,~] = trapveltraj(wayPoints,numSamples,'EndTime',T);
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
