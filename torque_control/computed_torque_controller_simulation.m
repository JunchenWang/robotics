
function computed_torque_controller_simulation

% robot = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
robot = convert_robot_tree(importrobot('E:\data\URDF\iiwa7\iiwa7.urdf'));
robot.TCP = [eye(3), [0, 0, 0.213]'; 0 0 0 1];
n = robot.dof;
u = udpport("byte");
ptp([-40, 70, 0, -80, 0, -60, 0]/180*pi);
% ptp([0, -60, 80, -100, -90, 0]/180*pi);

lineTo2(robot, [-500, 0,0], axang2rotm([0,1,0,pi/6]));

    function F = Wrench(t, y)
        F = zeros(6,7);
        if t < 3 && t > 1
%              F(:,end) = [0, 0, 0, 0, 0, -10]';
             F(:,4) = [0, 0, 0, 0, 0, 0]';
        end
    end

    function lineTo2(robot, t, R)
        if nargin < 3
            R = eye(3);
        end
        if isrow(t)
            t = t';
        end
 
        MassMatrix = @(t, y) [eye(n), zeros(n); zeros(n), mass_matrix(robot, y(1:n))];
        opts = odeset('Mass',MassMatrix,'OutputFcn',@odeplot);
        start = queryJoints;
        kesai = cal_kuka_kesai(start);
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
        Kd = diag([1, 1, 1, 1, 1, 1] * 5);
        Ki = diag([1, 1, 1, 1, 1, 1] * 2);
        Dn = diag(ones(1,7) * 2);
        Kn = diag(ones(1,7) * 4);
        y0 = [start, zeros(1,n)]';
        clear computed_torque_controller2;
        clear computed_torque_controller1;
        tao = @(t, y) computed_torque_controller2(robot, t,y, tforms,vel, acc, Kp, Ki, Kd, dt)...
                      + nullspace_impedance_controller(robot, t, y, tforms, kesai, Dn, Kn, dt);
        control_target = @(t, y) manipulator_dynamics_extForce(robot, tao, @Wrench,t, y); 
        [t,y] = ode45(control_target,[0, T + 1],y0,opts);
        plot(t, y(:,1:n)');
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

