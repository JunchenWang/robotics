function control_simulation
u = udpport("byte");
[mass, inertia, A, M, ME] = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
n = length(mass);
dynamics = @(t, y) manipulator_dynamics(mass, inertia, A, M, ME, @joint_torque, @ext_wrench, t, y);
MassMatrix = @(t, y) mass_matrix(mass, inertia, A, M, y(1:n));
opts = odeset('Mass',MassMatrix,'OutputFcn',@odeplot);
y0 = zeros(2*n,1);
y0(4) = pi / 3;

ptp(y0(1:n)');
% freq = 50;
% time = 10;
% tspan = linspace(0, time, freq * time);
tspan = [0, 10];
[t,y] = ode45(dynamics,tspan,y0,opts);
plot(t, y');

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
        if isempty(flag)
            setJoints(y(1:n));
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
        wayPoints = [start',jts'];
        Freq = 200;
        r = rateControl(Freq);
        T = max(abs(jts - start) / vel);
        numSamples = round(T * Freq) + 1;
        [q,qd,qdd,tSamples,pp] = trapveltraj(wayPoints,numSamples);
        for i = 1 : numSamples
            setJoints(q(:,i));
            waitfor(r);
        end
    end
    function joints = queryJoints
        % ";" 表示查询关节角
        writeline(u,"robot;","192.168.3.34",7755);
        s = readline(u);
        joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
        joints = mod(joints + pi, 2*pi) - pi;
    end

end
