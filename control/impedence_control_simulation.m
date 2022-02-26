
function impedence_control_simulation
global u;
global dt;
global kesai;
robot = read_dynamics_file('F:\robotics\urdf\iiwa7\dynamics.txt');
u = udpport("byte");
ptp([0, 70, 0, -80, 0, -60, 0]/180*pi);
start = queryJoints;
kesai = cal_kuka_kesai(start);
lineTo2(robot, [-400,0,0], axang2rotm([0,1,0,pi/2]));


function F = Wrench(t)
if t < 0.6 && t > 0.3
    F = [0, 0, 0, 0, 0, 0]';
else
    F = zeros(6,1);
end

function lineTo2(robot, t, R, vel)
if nargin < 4
    vel = 500;
end
if nargin < 3
    R = eye(3);
end
if isrow(t)
    t = t';
end

start = queryJoints;
pre_q = start;
q = pre_q;
kesai = cal_kuka_kesai(start);
Ts = forward_kin_kuka(start);
Ts(1:3,4) = Ts(1:3,4) / 1000;
Te = Ts*[R,t / 1000;0 0 0 1];
T = norm(t) / vel;
Freq = 500;
dt = 1 / Freq;
r = rateControl(Freq);
numSamples = round(T * Freq) + 1;
data = zeros(7, numSamples);
[s,sd,sdd] = trapveltraj([0, 1],numSamples, 'EndTime', T);
tSamples = linspace(0,T,numSamples);
[tforms,vel,~] = transformtraj(Ts,Te,[0 T],tSamples, 'TimeScaling', [s;sd;sdd]);
% plot(tSamples, reshape(tforms(1,4,:),[1, numSamples]));
M = cat(3, 0.01*eye(3), 0.01*eye(3));
B = cat(3, 2 * eye(3), 2 * eye(3));
K = cat(3, 50 * eye(3), 50 * eye(3));
for i = 1 : numSamples
    dist = 1e10;
    angles = [];
    qd = (q - pre_q) / dt;
    data(:,i) = qd';
    pre_q = q;
    y = [q,qd]';
    f = Wrench(tSamples(i));
    X = impedence_control(robot, tforms(:,:,i), vel(:,i), M, B, K, y, f, dt);
%     X = impedence_control(robot, Te, zeros(6,1), M, B, K, y, f, dt);
    for cfg1 = -1 : 2 : 1
        for cfg2 =  -1 : 2 : 1
            for cfg3 = -1 : 2 : 1
                cfg = [cfg1, cfg2, cfg3];
                ang = inverse_kin_kuka_kesai(X(1:3,1:3), X(1:3,4) * 1000, cfg, kesai, pre_q);
                if ~isempty(ang) && norm(ang-pre_q)<dist
                    dist = norm(ang-pre_q);
                    angles = ang;
                end
            end
        end
    end
    if isempty(angles) || limit_check_kuka(angles)
        plot(data');
        error('no solution');
    end
    setJoints(angles);
    q = angles;
    waitfor(r);
end
plot(data');
function ptp(jts, vel)
if nargin < 2
    vel = 1;
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

function joints = queryJoints
global u;
% ";" 表示查询关节角
writeline(u,"robot;","192.168.2.191",7755);
s = readline(u);
joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';

function gohome
ptp([0,0,0,0,0,0,0]);

function setJoints(jt)
global u;
cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
    ,jt(6), jt(7));
writeline(u,cmd,"192.168.2.191",7755);

