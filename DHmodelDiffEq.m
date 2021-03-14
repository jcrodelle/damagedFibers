function dxdt = DHmodelDiffEq(t,x,tVecA,avgFR_A,tVecC,avgFR_C, params)
%% these are the firing rates of three populations - passed in as vector x
fP = x(1);
fE = x(2);
fI = x(3);
gNMDA = x(4);

%% weights between populations - passed in as vector params 
gABP = params(1);  
gCP = params(2);
gEP = params(3);
gIP = params(4);
gCE = params(5);
gIE = params(6);
gABI = params(7);
%% response function parameters
alphaP = params(8);
betaP = params(9);
alphaE = params(10);
betaE = params(11);
alphaI = params(12);
betaI = params(13);
alphaG = params(14);
betaG = params(15);
maxG = params(16);

%% parameters for the firing rate differential equations
tauP = 0.001;
tauE = 0.01;
tauI = 0.02;
tauNMDA = 1.0; 

maxP = 50;
maxE = 60;
maxI = 80;

fA = interp1(tVecA, avgFR_A, t); % Interpolate the data set (ft, f) at times t
fC = interp1(tVecC, avgFR_C, t); % Interpolate the data set (gt, g) at times t

if t > 0.99 
    fA = avgFR_A(end);
    fC = avgFR_C(end);
end


%% response curves
Winf = @(c) maxP*0.5*(1 + tanh((c-betaP)/alphaP)); 
Einf = @(c) maxE*0.5*(1 + tanh((c-betaE)/alphaE)); 
Iinf = @(c) maxI*0.5*(1 + tanh((c-betaI)/alphaI))+1; 
Minf = @(c) maxG*0.5*(1 + tanh((c-betaG)/alphaG));

%% differential equations for firing rates
dfWdt = (Winf(gABP*fA + (gCP + gNMDA)*fC  + gEP*fE - gIP*fI) - fP)/tauP;
dfEdt = (Einf(gCE*fC  - gIE*fI)-fE)/tauE;
dfIdt = (Iinf(gABI*fA) - fI)/tauI;
dgNMDAdt = (maxG*Minf(fP) - gNMDA)/tauNMDA;

dxdt = [dfWdt;dfEdt;dfIdt;dgNMDAdt];