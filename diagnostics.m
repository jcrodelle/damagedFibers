% diagnostics to calculate variables for quantifying pain response
% input:
% t_star_vec - vector of times the average firing rate of W neurons crosses
% threshold, f_star
% A0 - area under fW between t0 and tfin of t_star_vec (A_total in paper)
% A_star - area above threshold, between t0 and tfin
% Q_star - average pain response (called pi* in paper) = A_star/(tfin - t0)
% fmax - maximum firing rate of W neurons
% numCrossing - number of times fW crosses threshold.

% output:
% t0_star - vector with avg and std for the FIRST time crosses thresh
% tN_star - vector with avg and std for the LAST time crosses thresh
% all of the others are also vectors with mean and std

% bigW_vec holds fW per realization fW_1 | fW_2 | ... 
function [t0_star, tN_star, A0, A_star, Q_star, fmax, numCrossings] = diagnostics(bigW_vec,f_star,t)

[RR,N] = size(bigW_vec);  %number of columns gives number of realizations

fmax_vec = zeros(RR,1);  %hold max fW for each realization
A0_vec = zeros(RR,1); % area under fW (entire thing)
A_star_vec = zeros(RR,1);  %area above f_star
Q_star_vec =  zeros(RR,1);  %avg pain response
t_star_bigVec = zeros(RR,50); %matrix to hold t_star_vec for each realization
numCrossings_vec = zeros(RR,1); % number of crossings for t_star_vec
t0_star_vec = zeros(RR,1); %first time fW crosses thresh
tN_star_vec = zeros(RR,1); %last time fW crosses thresh

% for each realization:
for i = 1:RR
    fW = bigW_vec(i,:);  
    [X,I] = max(fW); 
    fmax_vec(i) = X(1);
    % area under fW
    A0_vec(i) = trapz(t,fW);  % t in seconds
    
    %crossings:
    %subtract off the threshold:
    fW_aboveThresh = fW - f_star*ones(size(fW));
    B = fW_aboveThresh; % this is the subtracted firing rate
      % crossings - find all times this changes sign
    I = find(diff(B>=0));
    if isempty(I) %this did not cross threshold
        disp(['RR = ', num2str(i), ' No crossings!'])
        t0_star_vec(i) = 0;
        tN_star_vec(i) = 0;
        Q_star_vec(i) = 0;
        A_star_vec(i) = 0;
        numCrossings_vec(i) = 0;
    else
        t_temp = t(I);
        t_star_bigVec(i,1:length(I))=t(I);
        numCrossings_vec(i) = length(I);
        t0_star_vec(i) = t(I(1));
        tN_star_vec(i) = t(I(end));
        
        % area above 25:
        fW_aboveThresh(fW_aboveThresh<0) = 0;% set negative terms to 0
        
        A_star_vec(i) = trapz(t,fW_aboveThresh);
        
        % avg pain response:
        Q_star_vec(i) = A_star_vec(i)/(t_temp(end) - t_temp(1));
    end    
end  %end loop over realizations
% calculate avg and std for each
A0 = [mean(A0_vec) std(A0_vec)]; 
A_star = [mean(A_star_vec)  std(A_star_vec)];
Q_star = [mean(Q_star_vec) std(Q_star_vec)];
fmax = [mean(fmax_vec) std(fmax_vec)];
numCrossings = [mean(numCrossings_vec) std(numCrossings_vec)];
t0_star= [mean(t0_star_vec) std(t0_star_vec)];
tN_star = [mean(tN_star_vec) std(tN_star_vec)];

end
