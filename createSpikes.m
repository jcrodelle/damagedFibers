% will create and injure the spike trains for each fiber:
function [totalSpikeMat, totalTime, smoothM] = createSpikes(stimFiber, tFin)
baselineFreq = stimFiber(1);
stimFreq = stimFiber(2);
numberFibers = stimFiber(3);
initialBaseDuration = stimFiber(4);
durationOfStimulus= stimFiber(5);

%% baseline
% baseline for time prior to the fiber spiking: 
[spikeMat_baseline, tVec_baseline, spikeTimes_baseline] = poissonSpikeGen(baselineFreq, numberFibers,initialBaseDuration);
totalTime = tVec_baseline; % put this baseline into the total time vector
totalSpikeTimes = spikeTimes_baseline;  % put these baseline spike times into vector
totalSpikeMat = spikeMat_baseline;  % these are 0s and 1s

%% stimulus time
[spikeMat_stim, tVec_stim, spikeTimes_stim] = poissonSpikeGen(stimFreq, numberFibers, durationOfStimulus);
tVec_stim = (tVec_stim+initialBaseDuration);
totalTime = [totalTime tVec_stim];
spikeTimes_stim = spikeTimes_stim + initialBaseDuration;
totalSpikeTimes = [totalSpikeTimes spikeTimes_stim]; 
totalSpikeMat = [totalSpikeMat spikeMat_stim];
durationRemainderBaseline = tFin - (durationOfStimulus + initialBaseDuration); 

[spikeMat_remainderBaseline, tVec_remainderBaseline, spikeTimes_remainderBaseline] = poissonSpikeGen(baselineFreq, numberFibers,durationRemainderBaseline);
tVec_remainderBaseline = tVec_remainderBaseline + durationOfStimulus + initialBaseDuration;
spikeTimes_remainderBaseline = spikeTimes_remainderBaseline + durationOfStimulus + initialBaseDuration;
totalTime = [totalTime tVec_remainderBaseline];
totalSpikeTimes = [totalSpikeTimes spikeTimes_remainderBaseline];
totalSpikeMat = [totalSpikeMat spikeMat_remainderBaseline];

% change this into a smooth firing rate:
numbSpikesVec = sum(totalSpikeMat,1);
totalFR = (1000/(numberFibers)).*numbSpikesVec;
smoothWidth = 10;
smoothM = smooth(totalFR,smoothWidth,'rloess');
