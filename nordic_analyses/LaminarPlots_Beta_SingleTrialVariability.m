clc; clear all; close all;

subj       = 'S2_MN_NG_PC';     % Change participant number accordingly
datafolder = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\';

% outputfolder
folderOut = [datafolder,subj,filesep,'ResultsScriptVolume',filesep];
if ~exist(folderOut)
    mkdir(folderOut);
end

dirdata = ['D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\',subj,'\ResultsScriptVolume\'];
dirpoi = ['D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\',subj,'\ROI\'];

nrMaps = 880;

Beta_Std = xff([dirdata, subj, '_Betas_SingleTrial_Comp_H_bn_11-DepthLevels.smp']);
Beta_NOR = xff([dirdata, subj, '_Betas_SingleTrial_Comp_H_an_11-DepthLevels.smp']);
Beta_NN = xff([dirdata, subj, '_Betas_SingleTrial_Comp_H_an_NoNoise_11-DepthLevels.smp']);


poi_LH = xff([dirpoi, subj, '_ROI_LH_no_overlap.poi']);                    % Load the poi of one hemisphere.
HG_LH = poi_LH.POI(1).Vertices;   

Betas = zeros( size(HG_LH,1), 880);
Betas_NOR = zeros( size(HG_LH,1), 880);
Betas_NN = zeros( size(HG_LH,1), 880);

for i = 1:nrMaps
    temp = Beta_Std.Map(i).SMPData(HG_LH);
    Betas(:,i) = temp;
    temp2 = Beta_NOR.Map(i).SMPData(HG_LH);
    Betas_NOR(:,i) = temp2;
    temp3 = Beta_NN.Map(i).SMPData(HG_LH);
    Betas_NN(:,i) = temp3;
end 

D1 = linspace(1,870,80);
D2 = linspace(2,871,80);
D3 = linspace(3,872,80);
D4 = linspace(4,873,80);
D5 = linspace(5,874,80);
D6 = linspace(6,875,80);
D7 = linspace(7,876,80);
D8 = linspace(8,877,80);
D9 = linspace(9,878,80);
D10 = linspace(10,879,80);
D11 = linspace(11,880,80);

%% Depth 1
D1_Betas = Betas(:,D1);
D1_mean = mean(D1_Betas);
D1_se = std(D1_mean)/sqrt(length(D1_mean));

D1_Betas_NOR = Betas_NOR(:,D1);
D1_mean_NOR = mean(D1_Betas_NOR);
D1_se_NOR = std(D1_mean_NOR)/sqrt(length(D1_mean_NOR));

D1_Betas_NN = Betas_NN(:,D1);
D1_mean_NN = mean(D1_Betas_NN);
D1_se_NN = std(D1_mean_NN)/sqrt(length(D1_mean_NN));

%% Depth 2
D2_Betas = Betas(:,D2);
D2_mean = mean(D2_Betas);
D2_se = std(D2_mean)/sqrt(length(D2_mean));

D2_Betas_NOR = Betas_NOR(:,D2);
D2_mean_NOR = mean(D2_Betas_NOR);
D2_se_NOR = std(D2_mean_NOR)/sqrt(length(D2_mean_NOR));

D2_Betas_NN = Betas_NN(:,D2);
D2_mean_NN = mean(D2_Betas_NN);
D2_se_NN = std(D2_mean_NN)/sqrt(length(D2_mean_NN));

%% Depth 3
D3_Betas = Betas(:,D3);
D3_mean = mean(D3_Betas);
D3_se = std(D3_mean)/sqrt(length(D3_mean));

D3_Betas_NOR = Betas_NOR(:,D3);
D3_mean_NOR = mean(D3_Betas_NOR);
D3_se_NOR = std(D3_mean_NOR)/sqrt(length(D3_mean_NOR));

D3_Betas_NN = Betas_NN(:,D3);
D3_mean_NN = mean(D3_Betas_NN);
D3_se_NN = std(D3_mean_NN)/sqrt(length(D3_mean_NN));

%% Depth 4
D4_Betas = Betas(:,D4);
D4_mean = mean(D4_Betas);
D4_se = std(D4_mean)/sqrt(length(D4_mean));

D4_Betas_NOR = Betas_NOR(:,D4);
D4_mean_NOR = mean(D4_Betas_NOR);
D4_se_NOR = std(D4_mean_NOR)/sqrt(length(D4_mean_NOR));

D4_Betas_NN = Betas_NN(:,D4);
D4_mean_NN = mean(D4_Betas_NN);
D4_se_NN = std(D4_mean_NN)/sqrt(length(D4_mean_NN));

%% Depth 5
D5_Betas = Betas(:,D5);
D5_mean = mean(D5_Betas);
D5_se = std(D5_mean)/sqrt(length(D5_mean));

D5_Betas_NOR = Betas_NOR(:,D5);
D5_mean_NOR = mean(D5_Betas_NOR);
D5_se_NOR = std(D5_mean_NOR)/sqrt(length(D5_mean_NOR));

D5_Betas_NN = Betas_NN(:,D5);
D5_mean_NN = mean(D5_Betas_NN);
D5_se_NN = std(D5_mean_NN)/sqrt(length(D5_mean_NN));

%% Depth 6
D6_Betas = Betas(:,D6);
D6_mean = mean(D6_Betas);
D6_se = std(D6_mean)/sqrt(length(D6_mean));

D6_Betas_NOR = Betas_NOR(:,D6);
D6_mean_NOR = mean(D6_Betas_NOR);
D6_se_NOR = std(D6_mean_NOR)/sqrt(length(D6_mean_NOR));

D6_Betas_NN = Betas_NN(:,D6);
D6_mean_NN = mean(D6_Betas_NN);
D6_se_NN = std(D6_mean_NN)/sqrt(length(D6_mean_NN));

%% Depth 7
D7_Betas = Betas(:,D7);
D7_mean = mean(D7_Betas);
D7_se = std(D7_mean)/sqrt(length(D7_mean));

D7_Betas_NOR = Betas_NOR(:,D7);
D7_mean_NOR = mean(D7_Betas_NOR);
D7_se_NOR = std(D7_mean_NOR)/sqrt(length(D7_mean_NOR));

D7_Betas_NN = Betas_NN(:,D7);
D7_mean_NN = mean(D7_Betas_NN);
D7_se_NN = std(D7_mean_NN)/sqrt(length(D7_mean_NN));

%% Depth 8
D8_Betas = Betas(:,D8);
D8_mean = mean(D8_Betas);
D8_se = std(D8_mean)/sqrt(length(D8_mean));

D8_Betas_NOR = Betas_NOR(:,D8);
D8_mean_NOR = mean(D8_Betas_NOR);
D8_se_NOR = std(D8_mean_NOR)/sqrt(length(D8_mean_NOR));

D8_Betas_NN = Betas_NN(:,D8);
D8_mean_NN = mean(D8_Betas_NN);
D8_se_NN = std(D8_mean_NN)/sqrt(length(D8_mean_NN));

%% Depth 9
D9_Betas = Betas(:,D9);
D9_mean = mean(D9_Betas);
D9_se = std(D9_mean)/sqrt(length(D9_mean));

D9_Betas_NOR = Betas_NOR(:,D9);
D9_mean_NOR = mean(D9_Betas_NOR);
D9_se_NOR = std(D9_mean_NOR)/sqrt(length(D9_mean_NOR));

D9_Betas_NN = Betas_NN(:,D9);
D9_mean_NN = mean(D9_Betas_NN);
D9_se_NN = std(D9_mean_NN)/sqrt(length(D9_mean_NN));

%% Depth 10
D10_Betas = Betas(:,D10);
D10_mean = mean(D10_Betas);
D10_se = std(D10_mean)/sqrt(length(D10_mean));

D10_Betas_NOR = Betas_NOR(:,D10);
D10_mean_NOR = mean(D10_Betas_NOR);
D10_se_NOR = std(D10_mean_NOR)/sqrt(length(D10_mean_NOR));

D10_Betas_NN = Betas_NN(:,D10);
D10_mean_NN = mean(D10_Betas_NN);
D10_se_NN = std(D10_mean_NN)/sqrt(length(D10_mean_NN));

%% Depth 11
D11_Betas = Betas(:,D11);
D11_mean = mean(D11_Betas);
D11_se = std(D11_mean)/sqrt(length(D11_mean));

D11_Betas_NOR = Betas_NOR(:,D11);
D11_mean_NOR = mean(D11_Betas_NOR);
D11_se_NOR = std(D11_mean_NOR)/sqrt(length(D11_mean_NOR));

D11_Betas_NN = Betas_NN(:,D11);
D11_mean_NN = mean(D11_Betas_NN);
D11_se_NN = std(D11_mean_NN)/sqrt(length(D11_mean_NN));

%% Vector of mean responses across 11 depths
x = [mean(D1_mean) mean(D2_mean) mean(D3_mean) mean(D4_mean) mean(D5_mean) mean(D6_mean) mean(D7_mean) mean(D8_mean) mean(D9_mean) mean(D10_mean) mean(D11_mean)];
std_dev = [D1_se D2_se D3_se D4_se D5_se D6_se D7_se D8_se D9_se D10_se D11_se];

y = [mean(D1_mean_NOR) mean(D2_mean_NOR) mean(D3_mean_NOR) mean(D4_mean_NOR) mean(D5_mean_NOR) mean(D6_mean_NOR) mean(D7_mean_NOR) mean(D8_mean_NOR) mean(D9_mean_NOR) mean(D10_mean_NOR) mean(D11_mean_NOR)];
std_dev_NOR = [D1_se_NOR D2_se_NOR D3_se_NOR D4_se_NOR D5_se_NOR D6_se_NOR D7_se_NOR D8_se_NOR D9_se_NOR D10_se_NOR D11_se_NOR];

z = [mean(D1_mean_NN) mean(D2_mean_NN) mean(D3_mean_NN) mean(D4_mean_NN) mean(D5_mean_NN) mean(D6_mean_NN) mean(D7_mean_NN) mean(D8_mean_NN) mean(D9_mean_NN) mean(D10_mean_NN) mean(D11_mean_NN)];
std_dev_NN = [D1_se_NN D2_se_NN D3_se_NN D4_se_NN D5_se_NN D6_se_NN D7_se_NN D8_se_NN D9_se_NN D10_se_NN D11_se_NN];


figure() % Do we want to plot deep to superficial or superficial to deep?
hold on
delta = .05; % Adjust manually
h = plot(x);
set(h, 'Color',[0 0.4470 0.7410], 'LineWidth', 3)

p = plot(y);
set(p, 'Color',[0.8500 0.3250 0.0980], 'LineWidth', 3)

q = plot(z);
set(q, 'Color',[0.9290 0.6940 0.1250], 'LineWidth', 3)

errorbar((1:numel(x))-delta, x,std_dev,'o','linewidth',2,'MarkerSize',2', 'MarkerFaceColor', [0 0.4470 0.7410], 'Color',[0 0.4470 0.7410]) % Add X input
errorbar((1:numel(y)),y,std_dev_NOR,'o','linewidth',2,'MarkerSize',2', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'Color', [0.8500 0.3250 0.0980]) % Add X input
errorbar((1:numel(z))+delta,z,std_dev_NN,'o','linewidth',2,'MarkerSize',2', 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'Color', [0.9290 0.6940 0.1250]) % Add X input

legend('Ori', 'Nor\_def', 'Nor\_nn') % Location of legend somewhere else

f = gca;

f.XLim = ([1-delta 11+delta]);

h = gcf;
%saveas(h,[folderOut, subj, '_LaminarPlot.fig']);
%saveas(h, [folderOut, subj,'_LaminarPlot.svg']);
%saveas(h, [folderOut, subj, '_LaminarPlot.png']);
%save([folderOut, [subj, '_LaminarPlot.mat']],'x', 'y', 'z','std_dev', 'std_dev_NOR', 'std_dev_NN');


