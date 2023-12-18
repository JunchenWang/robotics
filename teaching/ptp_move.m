function ptp_move(port, jts, vel)
if nargin < 3
    vel = 2;
end
start = query_joints(port);
wayPoints = [start(:),jts(:)];
Freq = 200;
rate = rateControl(Freq);
T = max(abs(wayPoints(:,1) - wayPoints(:,2)) / vel);
numSamples = round(T * Freq) + 1;
jt = trapveltraj(wayPoints,numSamples);
for i = 1 : numSamples
    set_joints(port, jt(:,i));
    waitfor(rate);
end
end

