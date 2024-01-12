%%% Create figure 2

%% set path
datafolder = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\';
subj       = 'S2_MN_NG_PC';     % Change participant number accordingly

%% Load data for figures
SubP1 = load([datafolder, subj, '\ResultsScriptVolume\',subj,'_SpatialCorrelation_BetasAcrossCond.mat']);
SubP2 = load([datafolder, subj, '\ResultsScriptVolume\',subj,'_Run2runCorrelation_perMethod.mat']);
SubP3 = load([datafolder, subj, '\ResultsScriptVolume\',subj,'_SpatialCorrelation_BetasHL.mat']);

%% Subplot of Spatial Correlations of Betas for all conditions
CorrBetas = SubP1.CorrBetas;

groupCorrBetas = 0.*CorrBetas;  groupCorrBetas(1,:,:)  = 1; groupCorrBetas(2,:,:) = 2;
groupCorrBetas2 = 0.*CorrBetas; groupCorrBetas2(:,1,:) = 1; groupCorrBetas2(:,2,:) = 2; groupCorrBetas2(:,3,:) = 3; groupCorrBetas2(:,4,:) = 4;
groupCorrBetas2(:,5,:) = 5; groupCorrBetas2(:,6,:) = 6;

figure('Color',[1,1,1],'Position',[281 19 1160 500]),
subplot(1,14,1:6)
mb = boxchart(groupCorrBetas2(:),CorrBetas(:), 'GroupByColor', groupCorrBetas(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');
mb(1).BoxFaceColor = [144 103 167]./256;
mb(2).BoxFaceColor = [0.4660 0.6740 0.1880];
ylabel('Spatial Correlation', 'FontName', 'Arial', 'FontSize', 12)
title('Spatial Correlation per Condition', 'FontName', 'Arial', 'FontSize', 12)
%xlabel('Conditions','interpreter','Latex')
xticks([1, 2, 3, 4, 5, 6])
xticklabels({'PredH', 'MispredH', 'UnpredH', 'PredL', 'MispredL', 'UnpredL'})
legend({'NORdef vs Original', 'NORnn vs Original'}, 'FontName', 'Arial', 'FontSize', 10)
ylim([0.3 1])
xlim([0 7])

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';



%% Subplot of Pairwise Correlations

PPCor_BnH = SubP2.PPCor_BnH;
PPCor_AnH = SubP2.PPCor_AnH;
PPCor_AnNSH = SubP2.PPCor_AnNSH;

methCount = [repmat(1,size(PPCor_BnH,1),1); repmat(2,size(PPCor_BnH, 1),1); repmat(3,size(PPCor_BnH, 1),1)];

subplot(1,14,8:10)
boxchart([PPCor_BnH; PPCor_AnH; PPCor_AnNSH], 'GroupByColor', methCount, 'WhiskerLineStyle' ,'none','MarkerStyle','none')
title('Run to run Correlation', 'FontName', 'Arial', 'FontSize', 12)
xlabel('Dataset', 'FontName', 'Arial', 'FontSize', 12)
ylabel('Correlation', 'FontName', 'Arial', 'FontSize', 12)
set(gca, 'XTick', [])
%legend('Standard','NORDIC', 'NORDIC NN')
ylim([0.1 0.35])

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';


%% Subplot of Run to Average BN - LORO

CorrHL = SubP3.CorrHL;
CorrH = [CorrHL(1,:); CorrHL(3,:); CorrHL(5,:)];
meth = SubP3.meth;
meth = [meth(1,:); meth(3,:); meth(5,:)];

subplot(1,14,12:14)
mb = boxchart(CorrH(:),  'GroupByColor', meth(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');

ylabel('Cross Validated Correlation', 'FontName', 'Arial', 'FontSize', 12)
title('Correlation single run to Average', 'FontName', 'Arial', 'FontSize', 12)
%xticklabels(['PredH'; 'PredL'; 'PredH'; 'PredL'; 'PredH'; 'PredL'])
%xticks([0.66, 1.66, 3, 4, 5.33, 6.33])
legend({'Original', 'NORdef', 'NORnn'}, 'FontName', 'Arial', 'FontSize', 10)
ylim([0.25 0.5])
xlabel('Dataset', 'FontName', 'Arial', 'FontSize', 12)
set(gca, 'XTick', [])
ylim([0.25 0.45])


%xlim([0 7])

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';

%% Save the current figure

s = gcf;
folderOut = [datafolder, subj, '\ResultsScriptVolume\'];


saveas(s,[folderOut, [subj,'_Figure2Correlations.fig']]);
saveas(s, [folderOut, [subj,'_Figure2Correlations.svg']]);
%saveas(s, [folderOut, [subj,'_Figure2Correlations.png']]);

