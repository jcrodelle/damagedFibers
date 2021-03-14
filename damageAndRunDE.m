% this function will damage fibers and run the diff eq
function damageAndRunDE(numberRandomSeeds, Tfin, injury_param_A, injury_type_A, percentInjury_A,injury_param_C, injury_type_C, percentInjury_C, name)

%weights between populations:
gABP= 0.8;
gCP= 0.8;
gEP=0.35;
gIP=1.8;
gCE= 1.6; 
gIE = 0.6;
gABI = 0.8;

% response curves
alphaP= 11.5;
betaP = 28.2;
alphaE= 5.2;
betaE= 29.2;
alphaI= 9.5;
betaI= 28.0;
alphaG= 10.0;
betaG= 38.0;
maxG = 2.0;

% combine in one vector
params = [gABP,gCP,gEP,gIP,gCE,gIE,gABI,alphaP,betaP,alphaE,betaE, alphaI,betaI,alphaG,betaG,maxG];

% info for input fibers
baselineFreq = 1;     % baseline frequency on all fibers
numberFibers_A = 380; % number of AB fibers in bundle 
numberFibers_C = 820; % number of AD and C fibers in bundle
durationOfStimulus_A = 0.02; % AB response
durationOfBaseline_A = 0.5;  % baseline before pain
durationOfStimulus_C = 0.21; % C response
durationOfBaseline_C = 0.59; % baseline before pain
stimFreq_A = 40;  % - freq response to stimulus on AB fibers
stimFreq_C = 22;  % - freq response to stimulus on C fibers
    
% time vector
tspan = 0:0.001:Tfin;  
% information for each fiber:
stimInfoA = [baselineFreq, stimFreq_A, numberFibers_A, durationOfBaseline_A, durationOfStimulus_A];
stimInfoC = [baselineFreq, stimFreq_C, numberFibers_C, durationOfBaseline_C, durationOfStimulus_C];
% run over loops of random seeds for each firing rate frequency:
randSeed_vec = 1:numberRandomSeeds;

% will hold the firing rates for each realization
bigWvec = zeros(length(randSeed_vec), length(tspan));
bigEvec = zeros(length(randSeed_vec), length(tspan));
bigIvec = zeros(length(randSeed_vec), length(tspan));

bigWvec_normal = zeros(length(randSeed_vec), length(tspan));
bigEvec_normal = zeros(length(randSeed_vec), length(tspan));
bigIvec_normal = zeros(length(randSeed_vec), length(tspan));

for RR = 1:length(randSeed_vec)
    rng(RR);
    disp(['RR = ', num2str(RR)]);
    
    % create Poisson spike trains for NORMAL SIMULATION:
    [spikeMat_A, tVec_A, smoothFR_A] = createSpikes(stimInfoA, Tfin);
    [spikeMat_C, tVec_C, smoothFR_C] = createSpikes(stimInfoC, Tfin);
    
    
    % initital firing rates for all populations
    x0 = [0;0;0;0];
    
    % run ODE45 for NORMAL spike trains
    [t_normal,x_normal] = ode45(@(t,x) DHmodelDiffEq(t,x,tVec_A,smoothFR_A,tVec_C,smoothFR_C,params),tspan,x0);
    fW_normal = x_normal(:,1);
    fE_normal = x_normal(:,2);
    fI_normal = x_normal(:,3);
    
    % save NORMAL information into big 3d matrix
    bigWvec_normal(RR,:) = fW_normal;
    bigEvec_normal(RR,:) = fE_normal;
    bigIvec_normal(RR,:) = fI_normal;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INJURE Fibers:
    % injure A fibers
    if percentInjury_A > 0
        numberInjuries_A = round(percentInjury_A*numberFibers_A); %number of fibers to be injured
        I = randperm(numberFibers_A,numberInjuries_A); % index of fibers to be injured
        injuredSpikeMat_A = spikeMat_A(I,:); % injured fibers' spike vectors
        
        % introduce new matrix that contains the injured spike trains too
        newSpikeMat_C = spikeMat_A;
        for i = 1:numberInjuries_A
            foo = injuredSpikeMat_A(i,:); % injured spike vector
            injury_param = injury_param_A;
            injury_type = injury_type_A;
            [injuredTrain_A] = injureSpikeTrain(foo, injury_type, injury_param);
            newSpikeMat_A(I(i),:) = injuredTrain_A;
        end
        
        % after injury
        smoothWidth = 5.0; %ms
        numbSpikesVec_injured_A = sum(newSpikeMat_A,1);
        totalFR_injured_A = (1000/(numberFibers_A)).*numbSpikesVec_injured_A;
        smoothM_injured_A = smooth(totalFR_injured_A,smoothWidth,'rloess');
    else
        smoothM_injured_A = smoothFR_A;
    end
    
    if percentInjury_C>0
        numberInjuries_C = round(percentInjury_C*numberFibers_C); %number of fibers to be injured
        I = randperm(numberFibers_C,numberInjuries_C); % index of fibers to be injured
        injuredSpikeMat_C = spikeMat_C(I,:); % injured fibers' spike vectors
        
        % introduce new matrix that contains the injured spike trains too
        newSpikeMat_C = spikeMat_C;
        for i = 1:numberInjuries_C
            foo = injuredSpikeMat_C(i,:); % injured spike vector
            injury_type = injury_type_C;
            injury_param = injury_param_C;
            [injuredTrain_C] = injureSpikeTrain(foo, injury_type, injury_param);
  
            newSpikeMat_C(I(i),:) = injuredTrain_C;
%             size(newSpikeMat_C)
        end
        
        % after injury
        smoothWidth = 5.0; %ms
        numbSpikesVec_injured_C = sum(newSpikeMat_C,1);
        totalFR_injured_C = (1000/(numberFibers_C)).*numbSpikesVec_injured_C;
        smoothM_injured_C = smooth(totalFR_injured_C,smoothWidth);
    else
        smoothM_injured_C = smoothFR_C;
    end
    
    % fixing length of poisson output by adding last time point
    if length(tVec_C) ~= length(tspan)
        tVec_A = [tVec_A 1.0];
        tVec_C = [tVec_C 1.0];
        
        smoothFR_A = [smoothFR_A; smoothFR_A(end)];
        smoothFR_C = [smoothFR_C; smoothFR_C(end)];
        
        smoothM_injured_C = [smoothM_injured_C;  smoothM_injured_C(end)];
        smoothM_injured_A = [smoothM_injured_A;  smoothM_injured_A(end)];
    end
    % initital firing rates for populations
    x0 = [0;0;0;0];
    
    % run ODE45 for the diff eq in DHPPv1_MB
    [t,x_injure] = ode45(@(t,x) DHmodelDiffEq(t,x,tVec_A,smoothM_injured_A,tVec_C,smoothM_injured_C,params),tspan,x0);
    
    fW_injure = x_injure(:,1);
    fE_injure = x_injure(:,2);
    fI_injure = x_injure(:,3);
    gNMDA = x_injure(:,4);
    
    
    % save information into big 3d matrix
    bigWvec(RR,:) = fW_injure;
    bigEvec(RR,:) = fE_injure;
    bigIvec(RR,:) = fI_injure;
end
save([num2str(length(randSeed_vec)), 'realizations_',name,'_normal'],'bigWvec_normal','bigEvec_normal','bigIvec_normal')
save([num2str(length(randSeed_vec)), 'realizations_',name,'_injured'],'bigWvec','bigEvec','bigIvec')
end
