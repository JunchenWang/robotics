numSamples = 1000;
[q,qd,qdd,tSamples,pp] = trapveltraj([0, 1],numSamples);
plot(tSamples, q);