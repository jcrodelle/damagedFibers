% plot poplation responses
function plot_pain(numberRandomSeeds,name,percentInjury_C,datadir)

%avg normal, std normal, avg inj, std injure
Qstar_vec = zeros(1,4);
Astar_vec = zeros(1,4);
A0_vec = zeros(1,4);
fmax_vec = zeros(1,4);
numCrossings_vec = zeros(1,4);

% name = [num2str(numberRandomSeeds),'realizations',injury_type_C,'_percentInjure',num2str(100*percentInjury_C),endName];

newName = [num2str(numberRandomSeeds),'realizations_',name];
normalName = [newName,'_normal'];
injureName = [newName,'_injured'];

load([datadir, normalName],'bigWvec_normal','bigEvec_normal','bigIvec_normal')
load([datadir,injureName],'bigWvec','bigEvec','bigIvec')

% normal
avgW_normal = mean(bigWvec_normal,1); %row mean
avgE_normal = mean(bigEvec_normal,1);
avgI_normal = mean(bigIvec_normal,1);

stdW_normal = std(bigWvec_normal,0,1); %row mean
stdE_normal = std(bigEvec_normal,0,1);
stdI_normal = std(bigIvec_normal,0,1);

% injured:
avgW = mean(bigWvec,1); %row mean
avgE = mean(bigEvec,1);
avgI = mean(bigIvec,1);

stdW = std(bigWvec,0,1); %row mean
stdE = std(bigEvec,0,1);
stdI = std(bigIvec,0,1);

%%
tFin = 1;
dt = tFin/(length(avgW)-1);
t = 0:dt:tFin;


%% run diagnostics
f_star = 25; %threshold for firing rate to be painful
% normal neurons
[t0_star_normal, tN_star_normal, A0_normal, A_star_normal, Q_star_normal, fmax_normal, numCrossings_normal] = diagnostics(bigWvec_normal,f_star,t);

% injured neurons
[t0_star_injured, tN_star_injured, A0_injured, A_star_injured, Q_star_injured, fmax_injured, numCrossings_injured] = diagnostics(bigWvec,f_star,t);
%%
Qstar_vec(:) = [Q_star_normal Q_star_injured];
Astar_vec(:) = [A_star_normal A_star_injured];
A0_vec(:) = [A0_normal A0_injured];
fmax_vec(:) = [fmax_normal fmax_injured];
numCrossings_vec(:) = [numCrossings_normal numCrossings_injured];


figure
subplot(3,1,1)
shadedErrorBar(t,avgW_normal,stdW_normal,'lineprops','b')
hold on
shadedErrorBar(t,avgW,stdW,'lineprops','r')
hold on
plot(t,25*ones(size(t)),'k--')
axis([0 tFin 0 1.05*(max(avgW_normal)+max(stdW_normal))])
title([num2str(100*percentInjury_C),'injured neurons, Projection neurons (P)'])
set(gca,'FontSize',20.0)
subplot(3,1,2)
shadedErrorBar(t,avgE_normal,stdE_normal,'lineprops','g')
hold on
shadedErrorBar(t,avgE,stdE,'lineprops','r')
axis([0 tFin 0 1.05*(max(avgE_normal)+max(stdE_normal))])
title('Excitatory neurons (E)')
ylabel('Firing rate (Hz)')
set(gca,'FontSize',20.0)
subplot(3,1,3)
shadedErrorBar(t,avgI_normal,stdI_normal,'lineprops','k')
hold on
shadedErrorBar(t,avgI,stdI,'lineprops','r')
axis([0 tFin 0 1.05*(max(avgI_normal)+max(stdI_normal))])
title('Inhibitory neurons (I)')
xlabel('Time (s)')
set(gca,'FontSize',20.0)


G = figure;
shadedErrorBar(t,avgW_normal,stdW_normal,'lineprops','b')
hold on
shadedErrorBar(t,avgW,stdW,'lineprops','r')
hold on
plot(t,25*ones(size(t)),'k--')
axis([0 tFin 0 1.05*(max(avgW_normal)+max(stdW_normal))])
% title([num2str(100*percentInjury_C),'injured neurons, Projection neurons (W)'])
set(gca,'FontSize',25.0)
set(gcf, 'Position',  [300, 200, 500, 300])
saveas(G,['realizations_',name,'.png'])

%% 
F=figure;
errorbar(1,Qstar_vec(1), Qstar_vec(2),'sq','color',[0.08, 0.17, 0.55],'LineWidth',3.0,'MarkerFaceColor', 'b','MarkerSize',30.0)
hold on
errorbar(2, Qstar_vec(3), Qstar_vec(4),'sq','color',[0.64, 0.08, 0.18],'LineWidth',3.0,'MarkerFaceColor','r','MarkerSize',30.0)
title(['\pi^*'])
xlim([0.5 2.5])
set(gca,'xTick',[1, 2 ],'xTickLabel',{'Normal','Injured'},'FontSize',25.0)
set(gcf, 'Position',  [100, 100, 400, 300])
saveas(F,['Qstar_',name,'.png'])

F1 = figure;
errorbar(1,Astar_vec(1), Astar_vec(2),'sq','color',[0.08, 0.17, 0.55],'LineWidth',3.0,'MarkerFaceColor', 'b','MarkerSize',30.0)
hold on
errorbar(2, Astar_vec(3), Astar_vec(4),'sq','color',[0.64, 0.08, 0.18],'LineWidth',3.0,'MarkerFaceColor','r','MarkerSize',30.0)
title(['A^*'])
xlim([0.5 2.5])
set(gca,'xTick',[1, 2 ],'xTickLabel',{'Normal','Injured'},'FontSize',25.0)
set(gcf, 'Position',  [100, 100, 400, 300])
saveas(F1,['Astar_',name,'.png'])


F2 = figure;
errorbar(1,fmax_vec(1), fmax_vec(2),'sq','color',[0.08, 0.17, 0.55],'LineWidth',3.0,'MarkerFaceColor', 'b','MarkerSize',30.0)
hold on
errorbar(2, fmax_vec(3), fmax_vec(4),'sq','color',[0.64, 0.08, 0.18],'LineWidth',3.0,'MarkerFaceColor','r','MarkerSize',30.0)
title(['\pi_{max}'])
xlim([0.5 2.5])
set(gca,'xTick',[1, 2 ],'xTickLabel',{'Normal','Injured'},'FontSize',25.0)
set(gcf, 'Position',  [100, 100, 400, 300])
saveas(F2,['maxP_',name,'.png'])

end
