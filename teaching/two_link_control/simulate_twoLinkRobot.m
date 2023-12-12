function simulate_twoLinkRobot
u = udpport("byte");
robot = twoLinkRobot;
n = robot.dof;
pos0 = [pi/2, 0];
y0 = [pos0 + rand(1,2) * 0.001, 0, 0]';
tao = @(t, y) [0, 0]';
dynamic = @(t, y) twoLinkRobot_dynamics(t, y, robot, tao);
tspan = [0, 20];
opts = odeset('OutputFcn',@odeplot);
pret = 0;
[t, y] = ode45(dynamic, tspan, y0, opts);
plot(t, y(:,1:2)');
    function ret = odeplot(t, y, flag)
        if strcmp(flag, 'init') == 1
        elseif isempty(flag)
            Freq = 1 / (t(end) - pret);
            pret = t(end);
            r = rateControl(Freq);
            waitfor(r);
            setJoints(y(1:n, end));
        else
        end
        ret = 0;
    end
    function setJoints(jt)
        cmd = sprintf('robot;%f;%f;', jt(1), jt(2));
        writeline(u,cmd,"127.0.0.1",7755);
    end
end
