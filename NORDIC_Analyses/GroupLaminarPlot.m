%% This script creates the group laminar plot

%%% Variability is calculated as variance across participants %%%

dirdata = ['D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\ResultsGroupAnalysis\LaminarData\'];
subj = {'1', '2', '3', '5', '6', '7', '8', '9', '10', '11'};

for i = 1:size(subj,2)
    load([dirdata, 'S', subj{i}, '_MN_NG_PC_LaminarPlot.mat']);
    eval(['LaminarDataOri (i,1:11) = x']);
    eval(['LaminarDataNOR (i,1:11) = y']);
    eval(['LaminarDataNN (i,1:11) = z']);

end 

MeanOri = mean(LaminarDataOri);
MeanNOR = mean(LaminarDataNOR);
MeanNN = mean(LaminarDataNN);

SE_Ori = std(LaminarDataOri)/sqrt(length(LaminarDataOri));
SE_NOR = std(LaminarDataNOR)/sqrt(length(LaminarDataNOR));
SE_NN = std(LaminarDataNOR)/sqrt(length(LaminarDataNOR));



figure() 
hold on
delta = .05; % Adjust manually
h = plot(MeanOri);
set(h, 'Color',[0 0.4470 0.7410], 'LineWidth', 3)

p = plot(MeanNOR);
set(p, 'Color',[0.8500 0.3250 0.0980], 'LineWidth', 3)

q = plot(MeanNN);
set(q, 'Color',[0.9290 0.6940 0.1250], 'LineWidth', 3)

errorbar((1:numel(MeanOri))-delta, MeanOri,SE_Ori,'o','linewidth',2,'MarkerSize',2', 'MarkerFaceColor', [0 0.4470 0.7410], 'Color',[0 0.4470 0.7410]) % Add X input
errorbar((1:numel(MeanNOR)),MeanNOR,SE_NOR,'o','linewidth',2,'MarkerSize',2', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'Color', [0.8500 0.3250 0.0980]) % Add X input
errorbar((1:numel(MeanNN))+delta,MeanNN,SE_NN,'o','linewidth',2,'MarkerSize',2', 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'Color', [0.9290 0.6940 0.1250]) % Add X input

legend('Original', 'NORdef', 'NORnn', 'Location','northwest') % Location of legend somewhere else
hleg = legend;
hleg.FontSize = 10;
hleg.FontName = 'Arial';

f = gca;
f.XLim = ([1-delta 11+delta]);


xticks([1, 6, 11]);
xticklabels({'WM', 'GM', 'CSF'})

la = '$\textsf{$\beta$}$';
ylabel(la,'interpreter','Latex', 'fontsize', 12, 'fontname', 'Arial')

f.XAxis.FontSize = 12;
f.XAxis.FontName = 'Arial';
f.YAxis.FontSize = 12;
f.YAxis.FontName = 'Arial';

title('PredH', 'fontname', 'Arial', 'fontsize', 12)



h = gcf;
% 
saveas(h,[dirdata, 'GroupLaminarData_CompH.fig']);
saveas(h, [dirdata, 'GroupLaminarData_CompH.svg']);
saveas(h, [dirdata, 'GroupLaminarData_CompH.png']);

%% Check for Fede

meanDiff = MeanOri - MeanNOR;
varDiff = SE_Ori - SE_NOR;

figure
yyaxis left
plot(meanDiff)
hold on

yyaxis right
plot(varDiff)
legend('mean Difference', 'variance Difference')



