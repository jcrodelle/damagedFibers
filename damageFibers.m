%% Master file for damaging fibers
% you can damage either C or A fibers.
clear all; close all;

%% Set Basic Parameters
% Set number of random seeds (25 or 30)
numberRandomSeeds = 30;
% Set final time in seconds
Tfin = 1.0;

% make sure to change any parameters in the damage below
% flag for A fiber damge
flag_Adamage = 1; % 1 is yes to delay, 0 is no damage to A fibers
% flag for C fiber damage
% 0 is no damage,
% 1 is intermittant blocking,
% 2 is evoked potential,
% 3 is increased refractory,
flag_Cdamage = 3;

% set percentage injured
percentInjury_C = 0.25;
percentInjury_A = 0.5;

% Choose damage to A fibers (optional)
if flag_Adamage == 1
    injury_type_A = 'delay_train';
    injury_param_A.delay = 125.0;  %set b/t 100-300 ms
else
    percentInjury_A = 0.0; %set to 0 to damage only C fibers
    injury_param_A.delay = 0.0;
    injury_type_A = 'delay_train';
end

% name your data files:
beginName = [num2str(100*percentInjury_C),'percentInjuryC_',num2str(100*percentInjury_A),'percentInjuryA_'];


%% Choose to damage C fibers w/ intermittant blocking
if flag_Cdamage == 1
    injury_type_C = 'interm_blocking';
    % Set injury_params
    div = 20;  % pi/20
    injury_param_C.interm_freq =  pi/div; %100 ms
    % set the name of files
    endName = [injury_type_C,'_freqOver',num2str(div)];
    %increase refractory period
elseif flag_Cdamage == 2
    injury_type_C = 'evoke_potentials';
    %  Set injury_params
    injury_param_C.evokedProb = 0.3;   %prob of any spike evoking potentials
    injury_param_C.evoked = 2; %number of spikes generated by each spike
    injury_param_C.evokedSpace = 5; % space between generated spikes (ms)
    % set the name of files
    endName = [injury_type_C,'prob',num2str(100*injury_param_C.evokedProb), '_num',num2str(injury_param_C.evoked),'_space',num2str(injury_param_C.evokedSpace)];
elseif flag_Cdamage == 3
    injury_type_C = 'increase_refract';
    % Set injury_params
    injury_param_C.tau = 15;  %in ms
    % set the name of files
    endName = [injury_type_C,'tau',num2str(injury_param_C.tau)];
else
    percentInjury_C = 0.0;
    endName = 'noCinjury_Aonly';
end


%% Run the simulation
disp(['damaging ',num2str(100*percentInjury_C), '% of C fibers, and ', num2str(100*percentInjury_A), '% of A fibers']);
name = [beginName, endName];
damageAndRunDE(numberRandomSeeds, Tfin, injury_param_A, injury_type_A, percentInjury_A,injury_param_C, injury_type_C, percentInjury_C,name)


%% Plot the results from simulation
plot_pain(numberRandomSeeds, name, percentInjury_C,'')

