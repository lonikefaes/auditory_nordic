%%% Script to create Supplementary figure of Spatial correlations %%%

%%% It takes every individual Figure 2 and makes one big figure with all
%%% participants.

clc; clear all; close all;

%% Set main directories

datafolder = ['D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\'];

% outputfolder
folderOut = [datafolder,'ResultsGroupAnalysis',filesep];

allSubj = [1, 2, 3, 5, 6];
test = [1, 2, 3; 4, 5, 6; 7, 8, 9; 10, 11, 12; 13, 14, 15];%; 16, 17, 18; ...
    %19, 20, 21; 22, 23, 24; 25, 26, 27; 28, 29, 30];

figure('Color',[1,1,1],'Position',[281 19 1160 1200]),

%% set path
for i = 1:length(allSubj)
    subj       = ['S',num2str(allSubj(i)),'_MN_NG_PC'];     % Change participant number accordingly

    %% Load data for figures
    SubP1 = load([datafolder, subj, '\ResultsScriptVolume\',subj,'_SpatialCorrelation_BetasAcrossCond.mat']);
    SubP2 = load([datafolder, subj, '\ResultsScriptVolume\',subj,'_Run2runCorrelation_perMethod.mat']);
    SubP3 = load([datafolder, subj, '\ResultsScriptVolume\',subj,'_SpatialCorrelation_BetasHL.mat']);

    %% Subplot of Spatial Correlations of Betas for all conditions
    CorrBetas = SubP1.CorrBetas;

    groupCorrBetas = 0.*CorrBetas;  groupCorrBetas(1,:,:)  = 1; groupCorrBetas(2,:,:) = 2;
    groupCorrBetas2 = 0.*CorrBetas; groupCorrBetas2(:,1,:) = 1; groupCorrBetas2(:,2,:) = 2; groupCorrBetas2(:,3,:) = 3; groupCorrBetas2(:,4,:) = 4;
    groupCorrBetas2(:,5,:) = 5; groupCorrBetas2(:,6,:) = 6;

    subplot(5,3,test(i,1))
    hold on
    mb = boxchart(groupCorrBetas2(:),CorrBetas(:), 'GroupByColor', groupCorrBetas(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');
    mb(1).BoxFaceColor = [144 103 167]./256;
    mb(2).BoxFaceColor = [0.4660 0.6740 0.1880];
    % ylabel('Spatial Correlation','interpreter','Latex')
    %xlabel('Conditions','interpreter','Latex')
    xticks([1, 2, 3, 4, 5, 6])
    %xticklabels({'PredH', 'MispredH', 'UnpredH', 'PredL', 'MispredL', 'UnpredL'})
    %legend({'NORDIC vs Standard', 'NORDIC NN vs Standard' })
    ylim([0.3 1])


    b = gca;
    b.XAxis.FontSize = 8;
    b.XAxis.FontName = 'Arial';
    b.YAxis.FontSize = 8;
    b.YAxis.FontName = 'Arial';

    %% Subplot of Pairwise Correlations

    PPCor_BnH = SubP2.PPCor_BnH;
    PPCor_AnH = SubP2.PPCor_AnH;
    PPCor_AnNSH = SubP2.PPCor_AnNSH;

    methCount = [repmat(1,size(PPCor_BnH,1),1); repmat(2,size(PPCor_BnH, 1),1); repmat(3,size(PPCor_BnH, 1),1)];

    subplot(5,3,test(i,2))
    boxchart([PPCor_BnH; PPCor_AnH; PPCor_AnNSH], 'GroupByColor', methCount, 'WhiskerLineStyle' ,'none','MarkerStyle','none')
    %title('Run to run correlation for each method')
    xlabel('Dataset', 'FontName', 'Arial', 'FontSize', 12)
    %ylabel('Correlation')Comp
    set(gca, 'XTick', [])
    %legend('Standard','NORDIC', 'NORDIC NN')


    b = gca;
    b.XAxis.FontSize = 8;
    b.XAxis.FontName = 'Arial';
    b.YAxis.FontSize = 8;
    b.YAxis.FontName = 'Arial';

    %% Subplot of Run to Average BN - LORO

    CorrHL = SubP3.CorrHL;
    CorrH = [CorrHL(1,:); CorrHL(3,:); CorrHL(5,:)];
    meth = SubP3.meth;
    meth = [meth(1,:); meth(3,:); meth(5,:)];
    
    subplot(5,3,test(i,3))
    mb = boxchart(CorrH(:),  'GroupByColor', meth(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');

%     subplot(5,3,test(i,3))
%     mb = boxchart(methodID(:), CorrHL(:),  'GroupByColor', meth(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');

    %ylabel('Cross Validated Correlation','interpreter','Latex')
    %xticklabels(['PredH'; 'PredL'; 'PredH'; 'PredL'; 'PredH'; 'PredL'])
    %xticks([0.66, 1.66, 3, 4, 5.33, 6.33])
    %legend({'Standard', 'NORDIC', 'NORDIC NN'})
    ylim([0.0 0.5])
    set(gca, 'XTick', [])

    b = gca;
    b.XAxis.FontSize = 8;
    b.XAxis.FontName = 'Arial';
    b.YAxis.FontSize = 8;
    b.YAxis.FontName = 'Arial';

end

%% Save the current figure

s = gcf;
% folderOut = [datafolder, subj, '\ResultsScriptVolume\'];

saveas(s,[folderOut, 'SpatialCorr_S1-S5_Participants.fig']);
saveas(s, [folderOut, 'SpatialCorr_S1-S5_AllParticipants.svg']);
saveas(s, [folderOut, 'SpatialCorr_S1-S5_AllParticipants.png']);
