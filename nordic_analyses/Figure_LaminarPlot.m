%% Figure_LaminarPlot

clc; clear all; close all;

%% set path
datafolder = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\';

subjNR = [2, 7, 11];

folderOut = [datafolder,'ResultsGroupAnalysis',filesep];

figure('Color',[1,1,1],'Position',[281 19 1160 420]),

for i = 1:length(subjNR)
    subj       = ['S',num2str(subjNR(i)),'_MN_NG_PC'];     % Change participant number accordingly
    folderIn = [datafolder,subj(i),filesep,'ResultsScriptVolume',filesep];

    Lam = load([datafolder, subj, '\ResultsScriptVolume\',subj,'_LaminarPlot.mat']);

    delta = .05; % Adjust manually
    subplot(1,3,i)
    h = plot(Lam.x);
    set(h, 'Color',[0 0.4470 0.7410], 'LineWidth', 3)

    hold on

    p = plot(Lam.y);
    set(p, 'Color',[0.8500 0.3250 0.0980], 'LineWidth', 3)

    q = plot(Lam.z);
    set(q, 'Color',[0.9290 0.6940 0.1250], 'LineWidth', 3)
    

    errorbar((1:numel(Lam.x))-delta, Lam.x,Lam.std_dev,'o','linewidth',2,'MarkerSize',2', 'MarkerFaceColor', [0 0.4470 0.7410], 'Color',[0 0.4470 0.7410]) % Add X input
    errorbar((1:numel(Lam.y)),Lam.y,Lam.std_dev_NOR,'o','linewidth',2,'MarkerSize',2', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'Color', [0.8500 0.3250 0.0980]) % Add X input
    errorbar((1:numel(Lam.z))+delta,Lam.z,Lam.std_dev_NN,'o','linewidth',2,'MarkerSize',2', 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'Color', [0.9290 0.6940 0.1250]) % Add X input

    f = gca;

    f.XLim = ([1-delta 11+delta]);
    xticks([1, 6, 11]);
    xticklabels({'WM', 'GM', 'CSF'})

    ylabel('$\beta$  \% change','interpreter','Latex', 'fontsize', 16)

    title(['S', num2str(subjNR(i))], 'fontsize', 16)

end

legend('Original', 'NORdef', 'NORnn', 'Location', 'northwest') % Location of legend somewhere else

h = gcf; 

saveas(h,[folderOut, ['Figure11_LaminarPlots.fig']]);
saveas(h, [folderOut, ['Figure11_LaminarPlots.svg']]);
saveas(h, [folderOut, ['Figure11_LaminarPlots.png']]);

