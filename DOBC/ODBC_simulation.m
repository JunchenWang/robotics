function ODBC_simulation
u = udpport("byte");
robot = convert_robot_tree(importrobot('urdf\iiwa7\iiwa7.urdf'));
robot2 = robot;
robot2.mass = 1.2*robot.mass;% error
n = robot.dof;

kp = 25;
Kp = kp * eye(6);
Kd = 20 * eye(6);
Bn = diag(ones(1,7) * 4);
Kn = diag(ones(1,7) * 5);

MassMatrix = @(t, y) [eye(n), zeros(n, 2 * n); zeros(n), mass_matrix(robot, y(1:n)), zeros(n); zeros(n, 2 * n), eye(n)];
opts = odeset('Mass',MassMatrix,'OutputFcn',@odeplot);
y0 = zeros(3*n,1);
y0(1:n) = [0 75 0 -94 0 -81 0] / 180 * pi;
ptp(y0(1:n)');
kesai = cal_kuka_kesai(y0);
T0 = forward_kin_general(robot, y0);
refPose = T0;
% refPose(3,4) = refPose(3,4) + 0.5;
% refPose = trvec2tform([0.7 -.1 0.55]);
tInterval = [0, 10];
freq = 500;
N = tInterval(2) * freq + 1;
tSamples = linspace(tInterval(1), tInterval(2), N);
[scale,sd,sdd] = trapveltraj([0, 1], N, 'EndTime', tInterval(2));
[tforms,v,a] = transformtraj(refPose,refPose,tInterval,tSamples, 'TimeScaling', [scale;sd;sdd]);

Y = 10 * eye(n);
controller = @(t, y) DO_controller(robot2, tforms,v, a, Kp, Kd, Bn, Kn, kesai, freq, Y, t, y);      
control_target = @(t, y) manipulator_dynamics_observer(robot, controller, @Wrench, t, y); 
torque = [];
error = [];
speed = [];
tt = [];
[t,y] = ode15s(control_target,tInterval,y0,opts);
 X = forward_kin_general(robot, y(end,:));
        disp(logT(tform_inv(X) * refPose));
figure;
plot(t, y(:,2 * n + 1:end));
figure;
plot(tt, torque);
figure;
plot(tt, sqrt(sum(error(:,4:end).^2,2)));
figure;
plot(tt, speed);
    function F = Wrench(t, y)
        F = zeros(6,n);
        if t > 1
             F(:,4) = [0, 0, 0, 0, 0, 10]';
        end
    end

    function ret = odeplot(t, y, flag)
        if strcmp(flag, 'init') == 1
        elseif isempty(flag)
            setJoints(y(1:n, end));
            X = forward_kin_general(robot, y(1:n, end)) ; 
            Xd = tforms(:,:, round(t(end) * freq + 1));
            speed = [speed, y(n+1:2*n, end)];
            error = [error; twist_dist(X, Xd)];
            torque = [torque, controller(t(end), y(:,end))];
            tt = [tt, t(end)];
            % torque = [torque, y(2 * n+1:3*n, end)];
        else
        end
        ret = 0;
    end
    function setJoints(jt)
        cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
            ,jt(6), jt(7));
        writeline(u,cmd,"127.0.0.1",7755);
    end
    function ptp(jts, vel)
        if nargin < 2
            vel = .4;
        end
        start = queryJoints;
        wayPoints = [start',jts'];
        Freq = 200;
        r = rateControl(Freq);
        T = max(abs(jts - start) / vel);
        numSamples = round(T * Freq) + 1;
        jt = trapveltraj(wayPoints,numSamples);       
        for i = 1 : numSamples
            setJoints(jt(:,i));
            waitfor(r);
        end
    end
    function joints = queryJoints
        % ";" 表示查询关节角
        writeline(u,"robot;","127.0.0.1",7755);
        s = readline(u);
        joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
        joints = mod(joints + pi, 2*pi) - pi;
    end
end
