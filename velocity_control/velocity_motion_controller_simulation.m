
function velocity_motion_controller_simulation

% robot = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
robot = convert_robot_tree(importrobot('E:\data\URDF\iiwa7\iiwa7.urdf'));
n = robot.dof;
u = udpport("byte");
ptp([-40, 70, 0, -80, 0, -60, 0]/180*pi);
% ptp([0, -60, 80, -100, -90, 0]/180*pi);

lineTo2(robot, [-500, 0,0]); % axang2rotm([0,1,0,pi/2]));

    function lineTo2(robot, t, R)
        if nargin < 3
            R = eye(3);
        end
        if isrow(t)
            t = t';
        end
        start = queryJoints;
%         kesai = cal_kuka_kesai(start);
        Ts = forward_kin_general(robot, start);
        Te = Ts*[R,t / 1000;0 0 0 1];
        T = 5;
        Freq = 500;
        dt = 1 / Freq;
        numSamples = round(T * Freq) + 1;
        %data = zeros(robot.dof, numSamples);
        [s,sd,sdd] = trapveltraj([0, 1],numSamples, 'EndTime', T);
        tSamples = linspace(0,T,numSamples);
        [tforms,vel, acc] = transformtraj(Ts,Te,[0 T],tSamples, 'TimeScaling', [s;sd;sdd]);
        Kp = diag([100,100,100,100,100,100]*2);
        Kd = diag([1, 1, 1, 1, 1, 1] * 0); % PI controller
        Ki = diag([1, 1, 1, 1, 1, 1] * 2);
        q = start';
        qd = zeros(n, 1);
        clear motion_velocity_controller;
        for i = 1 : numSamples
            y = [q; qd];
            qd = motion_velocity_controller(robot, tSamples(i), y, tforms(:,:,i), vel(:,i), acc(:,i), Kp, Ki, Kd, dt);
            q = q + qd * dt;
            setJoints(q);
        end
        
    end
    function ret = odeplot(t, y, flag)
        if strcmp(flag, 'init') == 1
        elseif isempty(flag)
            setJoints(y(1:n, end));
%             torque = [torque, tao(t(end), y(:,end))];
        else
        end
        ret = 0;
    end
    function ptp(jts, vel)
        if nargin < 2
            vel = 1;
        end
        start = queryJoints;
        wayPoints = [start',jts'];
        Freq = 500;
        r = rateControl(Freq);
        T = max(abs(jts - start) / vel);
        numSamples = round(T * Freq) + 1;
        [q,~,~,~,~] = trapveltraj(wayPoints,numSamples);

        for i = 1 : numSamples
            setJoints(q(:,i));
            waitfor(r);
        end
    end

    function joints = queryJoints
        % ";" 表示查询关节角
        writeline(u,"robot;","127.0.0.1",7755);
        s = readline(u);
        joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
    end

    function gohome
        ptp([0,0,0,0,0,0,0]);
    end

    function setJoints(jt)
        cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
            ,jt(6), jt(7));
        writeline(u,cmd,"127.0.0.1",7755);
    end
end

