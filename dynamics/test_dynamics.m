function test_dynamics
[mass, inertia, A, M, ME] = read_dynamics_file('F:\MICR\MICSys\dynamics.txt');
u = udpport("byte");
robot = importrobot('E:\data\URDF\iiwa7\iiwa7.urdf');
robot.Gravity = [0, 0, -9.8];
robot.DataFormat = 'row';
gohome;
ptp([0, 70, 0, -80, 0, -60, 0]/180*pi);

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
tao = zeros(7, numSamples);
for i = 1 : numSamples
    setJoints(q(:,i));
    tao(:,i) = inverse_dynamics(mass, inertia, A, M, ME, q(:,i), qd(:,i), qdd(:,i), zeros(6,1));
%     tao(:,i) = inverseDynamics(robot, q(:,i)', qd(:,i)', qdd(:,i)');
    plot(1:i, tao(:,1:i)');
    waitfor(r);
end
% plot(1:i, qdd');
end
function joints = queryJoints
% ";" 表示查询关节角
writeline(u,"robot;","192.168.3.34",7755);
s = readline(u);
joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
end
function gohome
ptp([0,0,0,0,0,0,0]);
end
function setJoints(jt)
cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
    ,jt(6), jt(7));
writeline(u,cmd,"192.168.3.34",7755);
end
end

