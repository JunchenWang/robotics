
function computed_torque_setpoint_simulation

% robot = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
robot = convert_robot_tree(importrobot('E:\data\URDF\iiwa7\iiwa7.urdf'));
robot.TCP = [eye(3), [0, 0, 0.213]'; 0 0 0 1];
robot2 = robot;
robot2.mass(7) = robot2.mass(7) + 0.1;
n = robot.dof;
u = udpport("byte");
ptp([-40, 70, 0, -80, 0, -60, 0]/180*pi);
setpoint(robot, [-500, 0,0], axang2rotm([0,1,0,pi/6]));

    function F = Wrench(t, y)
        F = zeros(6,7);
        if t < 0.6 && t > 0.5
%              F(:,end) = [0, 0, 0, 0, 0, 10]';
%              F(:,4) = [0, 0, 0, 0, 0, 20]';
        end
    end

    function setpoint(robot, t, R)
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
        K = 10;
        Kp = diag([1,1,1,1,1,1]*K);
        Kd = diag([1, 1, 1, 1, 1, 1] * 2 * sqrt(K));
%         Ki = diag([1, 1, 1, 1, 1, 1] * 0);
        Dn = diag(ones(1,7) * 2 * sqrt(10));
        Kn = diag(ones(1,7) * 10);
        y0 = [start, zeros(1,n)]';
        tao = @(t, y) setpoint_controller(robot2, t,y, Te,zeros(6,1), zeros(6,1), Kp, Kd)...
                      + nullspace_controller(robot2, t, y, Te, kesai, Dn, Kn);
        control_target = @(t, y) manipulator_dynamics_fext(robot, tao, @Wrench,t, y); 
        [t,y] = ode45(control_target,[0, T],y0,opts);
        plot(t, y(:,1:n)');
        X = forward_kin_general(robot, y(end,:));
        disp(logT(tform_inv(X) * Te));
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

