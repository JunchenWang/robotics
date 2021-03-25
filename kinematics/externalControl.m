function externalControl
global u;
u = udpport("byte");
gohome;
numSamples = 501;
home = [0,0,0,0,0,0,0];
target = [-60,60,0,-60,0,-45,0];
wayPoints = [home',target'];
r = rateControl(200);
[q,qd,qdd,tSamples,pp] = trapveltraj(wayPoints,numSamples);
for i = 1 : numSamples
    setJoints(q(:,i));
    waitfor(r);
end


function gohome
global u;
writeline(u,"0;0;0;0;0;0;0;","127.0.0.1",7755);

function setJoints(jt)
global u;
cmd = sprintf('%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
        ,jt(6), jt(7));
writeline(u,cmd,"127.0.0.1",7755);  
