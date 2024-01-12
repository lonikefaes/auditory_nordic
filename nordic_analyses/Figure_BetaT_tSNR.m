%% Create figure 6

clc; clear all; close all;

%% set path
datafolder = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\';
subj       = 'S2_MN_NG_PC';     % Change participant number accordingly

folderOut = [datafolder,subj,filesep,'ResultsScriptVolume',filesep];

%% Load data for figures
tSNR = load([datafolder, subj, '\ResultsScriptVolume\',subj,'_tSNR.mat']);

figure('Color',[1,1,1],'Position',[281 19 500 420]),

it = 1;
temp          = squeeze(tSNR.BetasCond(it,1,:)) - squeeze(tSNR.BetasCond(it,2,:));
[Z bins]      = hist3([mean(tSNR.tsnrBN)',temp],[100,100]);
[xgrid ygrid] = meshgrid(bins{2}, bins{1} );
myp = plot(squeeze(tSNR.BetasCond(it,1,:)) - squeeze(tSNR.BetasCond(it,2,:)) , mean(tSNR.tsnrBN ) ,'.'); hold all
myp.Color = [0.301 0.75 0.933];
% [~, myc] = contour(xgrid,ygrid,Z,3, 'k-','LineWidth',2);xlim([-10 10]);axis square
% myc.LineColor = [0.3 0.3 0.3];
[my binsx ] = window_mean(  squeeze(tSNR.BetasCond(it,1,:)) - squeeze(tSNR.BetasCond(it,2,:)),  mean(tSNR.tsnrBN), linspace(7,60,50) );
plot(my,binsx,'r-','LineWidth',1)
xline(0)
xlabel('Betas Original - NORdef', 'FontName', 'Arial', 'FontSize', 12);
axis square;
ylabel('$\frac{Mean}{Std}$', 'interpreter','Latex','fontsize', 20)
xlim([-10 10]);

% b = gca;
% b.XAxis.FontSize = 12;
% b.XAxis.FontName = 'Arial';
% b.YAxis.FontSize = 12;
% b.YAxis.FontName = 'Arial';

h = gcf;

saveas(h,[folderOut, [subj,'_Figure6_BT_tSNR_NORdef.fig']]);
saveas(h, [folderOut, [subj,'_Figure6_BT_tSNR_NORdef.svg']]);


figure('Color',[1,1,1],'Position',[281 19 500 420]),


it = 1;
temp          = squeeze(tSNR.BetasCond(it,1,:)) - squeeze(tSNR.BetasCond(it,3,:));
[Z bins]      = hist3([mean(tSNR.tsnrBN)',temp],[100,100]);
[xgrid ygrid] = meshgrid(bins{2}, bins{1} );
myp = plot(squeeze(tSNR.BetasCond(it,1,:)) - squeeze(tSNR.BetasCond(it,3,:)) , mean(tSNR.tsnrBN ) ,'.'); hold all
myp.Color = [0.301 0.75 0.933];
% [~, myc] = contour(xgrid,ygrid,Z,3, 'k-','LineWidth',2);xlim([-10 10]);axis square
% myc.LineColor = [0.3 0.3 0.3];
[my binsx ] = window_mean(  squeeze(tSNR.BetasCond(it,1,:)) - squeeze(tSNR.BetasCond(it,3,:)),  mean(tSNR.tsnrBN), linspace(7,60,50) );
plot(my,binsx,'r-','LineWidth',1)
xline(0)
xlabel('Betas Original - NORnn', 'FontName', 'Arial', 'FontSize', 12);
axis square;
ylabel('$\frac{Mean}{Std}$', 'interpreter','Latex','fontsize', 20)
xlim([-10 10]);



% b = gca;
% b.XAxis.FontSize = 12;
% b.XAxis.FontName = 'Arial';
% b.YAxis.FontSize = 12;
% b.YAxis.FontName = 'Arial';


h = gcf;

saveas(h,[folderOut, [subj,'_Figure6_BT_tSNR_NORnn.fig']]);
saveas(h, [folderOut, [subj,'_Figure6_BT_tSNR_NORnn.svg']]);
