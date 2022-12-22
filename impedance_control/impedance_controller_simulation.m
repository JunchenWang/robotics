
function impedance_controller_simulation

robot = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
n = robot.dof;
u = udpport("byte");
ptp([-40, 70, 0, -80, 0, -60, 0]/180*pi);
% ptp([0, -60, 80, -100, -90, 0]/180*pi);

lineTo2(robot, [-500, 0,0]); % axang2rotm([0,1,0,pi/2]));

    function F = Wrench(t, y)
        if t < 3 && t > 1
            F = [10, 0, 0, 0, 0, 10]';
        else
            F = zeros(6,1);
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
        Kd = diag([100,100,100,1000,1000,1000]*0.1);
        Dd = diag([1, 1, 1, 10, 10, 10] * 5);
        y0 = [start, zeros(1,n)]';
        tao = @(t, y) impedance_controller2(robot, t,y, tforms,vel, acc, Dd, Kd, dt);
        control_target = @(t, y) manipulator_dynamics(robot, tao, @Wrench,t, y); 
        [t,y] = ode23(control_target,[0, T],y0,opts);
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

