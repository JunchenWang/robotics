function externalControl
global u;
u = udpport("byte");
ptp([0,60,0,-92,0,-60,0]/180*pi);
% lineTo3([-450,0,0]);
% lineTo3([0,0,-450]);
% lineTo3([450,0,0]);
% lineTo3([0,-450,0]);
for i = 1 : 1
lineTo([-450,0,0]);
lineTo([0,-450,0]);
lineTo([450,0,0]);
lineTo([0,450,0]);
end
 

function lineTo2(t, R, vel)
if nargin < 3
    vel = 200;
end
if nargin < 2
    R = eye(3);
end
if isrow(t)
    t = t';
end

start = queryJoints;
kesai = cal_kuka_kesai(start);
cfg=[sign(start(2)),sign(start(4)),sign(start(6))];
Ts = forward_kin_kuka(start);
Te = Ts*[R,t;0 0 0 1];
T = norm(t) / vel;
Freq = 200;
r = rateControl(Freq);
numSamples = round(T * Freq) + 1;
[s,sd,sdd] = trapveltraj([0, 1],numSamples, 'EndTime', T);
tSamples = linspace(0,T,numSamples);
[tforms,vel,~] = transformtraj(Ts,Te,[0 T],tSamples, 'TimeScaling', [s;sd;sdd]);
% plot(tSamples, reshape(tforms(1,4,:),[1, numSamples]));
oldang = [];
for i = 1 : numSamples
    oldkesai = kesai;
    [~, bd, ~] = inverse_kin_kuka(tforms(1:3,1:3, i), tforms(1:3,4, i), cfg);
%     if isempty(bd)
%         error('no solution');
%     end
%     [flag, a, b] = isInRange(bd, oldkesai);
%         if flag == 1
%             kesai = adjust_kesai(a, b, oldkesai);
%             angles = inverse_kin_kuka_kesai(tforms(1:3,1:3, i), tforms(1:3,4, i), cfg, kesai);
%         end
%      disp(kesai);
    flag = 0;
    if ~isempty(bd)
        [flag, a, b] = isInRange(bd, oldkesai);
        if flag == 1
            kesai = adjust_kesai(a, b, oldkesai);
            angles = inverse_kin_kuka_kesai(tforms(1:3,1:3, i), tforms(1:3,4, i), cfg, kesai);
        end
    end
    if flag == 0 || (~isempty(oldang) && norm(angles-oldang)>1)
        for cfg1 = -1 : 2 : 1
            for cfg2 =  -1 : 2 : 1
                for cfg3 = -1 : 2 : 1
                    cfg = [cfg1, cfg2, cfg3];
                    [~, bd, ~] = inverse_kin_kuka(tforms(1:3,1:3, i), tforms(1:3,4, i), cfg);
                    [flag, a, b] = isInRange(bd, oldkesai);
                    if flag == 1
                        kesai = adjust_kesai(a, b, oldkesai);
                        angles = inverse_kin_kuka_kesai(tforms(1:3,1:3, i), tforms(1:3,4, i), cfg, kesai);
                        if ~isempty(oldang) && norm(angles-oldang)<1
                            break;
                        else
                            flag = 0;
                        end
                    end
                end
                if flag == 1
                    break;
                end
            end
            if flag == 1
                break;
            end
        end
    end
    if flag == 0
        error('no solution');
    end
    if limit_check_kuka(angles)
        error('no solution');
    end
    disp(kesai);
    setJoints(angles);
    oldang = angles;
    waitfor(r);
end



function lineTo(t, R, vel)
fid = fopen('data.txt', 'a');
if nargin < 3
    vel = 100;
end
if nargin < 2
    R = eye(3);
end
if isrow(t)
    t = t';
end
d1 = 340;
d3 = 400;
d5 = 400;
d7 = 126;
start = queryJoints;
kesai = cal_kuka_kesai(start);
% kesai = cal_kuka_kesai(start);
% cfg=[sign(start(2)),sign(start(4)),sign(start(6))];
Ts = forward_kin_kuka(start);
Te = Ts*[R,t;0 0 0 1];
T = norm(t) / vel;
Freq = 200;
r = rateControl(Freq);
numSamples = round(T * Freq) + 1;
[s,sd,sdd] = trapveltraj([0, 1],numSamples, 'EndTime', T);
tSamples = linspace(0,T,numSamples);
[tforms,vel,~] = transformtraj(Ts,Te,[0 T],tSamples, 'TimeScaling', [s;sd;sdd]);
S = [0, 0, 1, 0, 0, 0;
    0, 1, 0, -d1, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, -1, 0, d1 + d3, 0, 0;
    0, 0, 1, 0, 0, 0;
     0, 1, 0, -(d1+d3+d5), 0, 0;
     0, 0, 1, 0, 0, 0];
% plot(tSamples, reshape(tforms(1,4,:),[1, numSamples]));
curJt = start;
step = tSamples(2) - tSamples(1);
for i = 1 : numSamples
    w = vel(1:3,i);
    v = vel(4:6, i);
    t = tforms(1:3,4, i);
    v = cross(w,-t) + v;
    Js = jacobian(S, curJt, 's');
    if cond(Js*Js') > 1e9
        error('near singularity');
    end
    vt = pinv(Js)*[w;v]*step; 
    curJt = curJt + vt';
    cfg = getConfig(curJt);
    curJt = inverse_kin_kuka_kesai(tforms(1:3,1:3, i), tforms(1:3,4, i), cfg, kesai);
    fprintf(fid, '%f %f %f %f %f %f %f\n', curJt);
    if limit_check_kuka(curJt)
        error('no solution');
    end
    setJoints(curJt);
    waitfor(r);
end
fclose(fid);

function lineTo3(t, R, vel)
if nargin < 3
    vel = 200;
end
if nargin < 2
    R = eye(3);
end
if isrow(t)
    t = t';
end
d1 = 340;
d3 = 400;
d5 = 400;
d7 = 126;
start = queryJoints;
theta3 = start(3);
R3 = [cos(theta3), -sin(theta3), 0;sin(theta3), cos(theta3), 0; 0, 0, 1];
Ts = forward_kin_kuka(start);
Te = Ts*[R,t;0 0 0 1];
T = norm(t) / vel;
Freq = 400;
r = rateControl(Freq);
numSamples = round(T * Freq) + 1;
[s,sd,sdd] = trapveltraj([0, 1],numSamples, 'EndTime', T);
tSamples = linspace(0,T,numSamples);
[tforms,vel,~] = transformtraj(Ts,Te,[0 T],tSamples, 'TimeScaling', [s;sd;sdd]);
S = [0, 0, 1, 0, 0, 0;
    0, 1, 0, -d1, 0, 0;
    [0, -1, 0] * R3', [d1 + d3, 0, 0]*R3';
    0, 0, 1, 0, 0, 0;
     [0, 1, 0] * R3', [-(d1+d3+d5), 0, 0] * R3';
     0, 0, 1, 0, 0, 0];
% plot(tSamples, reshape(tforms(1,4,:),[1, numSamples]));
curJt = [start(1), start(2), start(4:7)];
step = tSamples(2) - tSamples(1);
for i = 1 : numSamples
    w = vel(1:3,i);
    v = vel(4:6, i);
    t = tforms(1:3,4, i);
    v = cross(w,-t) + v;
    Js = jacobian(S, curJt, 's');
%     disp(cond(Js*Js'));
    if cond(Js) > 1e8
        error('near singularity');
    end
    vt = Js\[w;v]*step; 
    curJt = curJt + vt';
    cmdJt = [curJt(1:2), theta3, curJt(3:6)];
    kesai = cal_kuka_kesai(cmdJt);
    cfg = getConfig(cmdJt);
    cmdJt = inverse_kin_kuka_kesai(tforms(1:3,1:3, i), tforms(1:3,4, i), cfg, kesai);
    if limit_check_kuka(cmdJt)
        error('no solution');
    end
    setJoints(cmdJt);
    waitfor(r);
end


function moveTo(abc, t, cfg, vel)
if nargin < 4
    vel = .4;
end
[angles, ~, ~] = inverse_kin_kuka(EulZYX2R(abc), t, cfg);
if ~isempty(angles)
    ptp(angles, vel);
end

function ptp(jts, vel)
if nargin < 2
    vel = .4;
end
start = queryJoints;
wayPoints = [start',jts'];
Freq = 200;
r = rateControl(Freq);
T = max((jts - start) / vel);
numSamples = round(T * Freq) + 1;
[q,qd,qdd,tSamples,pp] = trapveltraj(wayPoints,numSamples);
for i = 1 : numSamples
    setJoints(q(:,i));
    waitfor(r);
end

function joints = queryJoints
global u;
writeline(u,"query;","127.0.0.1",7755);
s = readline(u);
joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f')';

function gohome
ptp([0,0,0,0,0,0,0]);

function setJoints(jt)
global u;
cmd = sprintf('%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
    ,jt(6), jt(7));
writeline(u,cmd,"127.0.0.1",7755);
