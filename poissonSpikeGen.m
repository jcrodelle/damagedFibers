function [spikeMat, tVec, spTimes] = poissonSpikeGen(fr, nTrials,t_end)
dt = 0.001; % 1 ms
nBins = round(t_end/dt);
spikeMat = rand(nTrials, nBins) < fr*dt;
tVec = 0:dt:t_end-dt;
spTimes = spikeMat.*tVec; % will contain zeros!



