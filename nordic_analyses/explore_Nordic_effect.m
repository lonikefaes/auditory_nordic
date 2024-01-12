clc;clear all;close all
%% set path
%addpath(genpath('C:\Program Files\NeuroElf_v10_5153\')); %neuroelf path
% addpath(genpath('D:\nordic\CircStat2012a'));% circular statistics
datafolder = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\';
% datafolder = 'D:\nordic\ToTin\ToTin\';

SaveBetaMapsMode = 1;    % 1 saves avergae betas per run per condition - 2 saves single trials (separated per condition in different maps)

subj       = 'S2_MN_NG_PC';

% outputfolder
folderOut = [datafolder,subj,filesep,'ResultsScriptVolume',filesep];
if ~exist(folderOut)
    mkdir(folderOut);
end

nruns      = 8;

xff(0,'transiosize','vtc',1048576);% read vtc fast
designCheck = 0;

colorsAll = [0 0 0.5; 0.3010 0.7450 0.9330; 0.4660, 0.6740, 0.1880];
colorsNOR = [0.3010 0.7450 0.9330; 0.4660, 0.6740, 0.1880];

%% Load ROIs from these names: function extract_roi_indices
roiAll    = {'TemporalLobe.msk'};
roinames    = {'HG_GM.msk','PP_GM.msk','PT_GM.msk','aSTG_GM.msk','pSTG_GM.msk'};
%        roinames    = {'HG_GM.msk'};
[roi_large, ~] = extract_roi_indices(roiAll,datafolder,subj);
[roi_dup, roiID] = extract_roi_indices(roinames,datafolder,subj);

%% Remove duplicates from roi's 
ii = roi_dup == roi_dup'; 
%numDup = sum(ii);
%Dup = find(numDup > 1);         % Find all the duplicates

[roi, IA, IC] = unique(roi_dup, 'first');
roiID = roiID(IA);

%% SLOW Load time series data: after and before nordic: reading time series
tsbn     = []; tsan = []; tsanNS = [];
Xbdiag   = [];  colnames = {};
runIdx = [];
locBN    = ~isnan(roi_large);locAN    = ~isnan(roi_large);  locANns = ~isnan(roi_large);
trial_x_run = [];
for itr    = 1:nruns
    disp(['Reading run ' , num2str(itr)]);
    % Load VTC before nordic of each run: function read_vtc which also
    % does  % signal change standardization
    [tsbn_t ,VTCInfo,tsnrBN(itr,:)] = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc' ],roi_large);
    locBN                           = locBN & sum(isnan(tsbn_t ))==0; % locate NaN in the data
    % load after nordic
    [tsan_t,~,  tsnrAN(itr,:)]      = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc' ],roi_large);
    locAN             = locAN & sum(isnan(tsan_t ))==0;% locate NaN in the data
    % load after nordic no noise
    [tsanNS_t,~,tsnrANns(itr,:)]   = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc' ],roi_large);
    locANns                          = locANns & sum(isnan(tsanNS_t ))==0;% locate NaN in the data
    % accumulate time series of each run in a big matrix
    tsbn      = [ tsbn ;   tsbn_t  ]; % make design matrix to have the right size
    tsan      = [ tsan ;   tsan_t  ]; % here we accumulate the data of all runs already standardized
    tsanNS    = [ tsanNS; tsanNS_t ];
    % creating design matrix here from prt
    X          = xff([datafolder,subj,filesep,'Prt_files',filesep,subj,'_run',num2str(itr),'_Cut_nogap_SingleTrialGLM.prt']);
    [X  , prednames]         = design_matrix_BV(X,VTCInfo.NrOfVolumes, VTCInfo.TR);% inside this function all the fields are defined to make X as BV
    colnames   = [colnames ,prednames]; % accumulate conditions x trials names accros runs

    Xbdiag     = blkdiag(Xbdiag,X);     % combine design matrix of each run in block diagonal structure
    ntimesrun(itr)  = size(X,1);
    trial_x_run = [ trial_x_run itr* ones(1, size(X,2) ) ];  % this vector defines the run for everytrial for be used in split half
end
%% nan excluding: we take only non nan voxels in the three datasets
tsbn   = tsbn(:,locANns & locAN & locBN);  % excluding NaN
tsan   = tsan(:,locANns & locAN & locBN);  % excluding NaN
tsanNS = tsanNS(:,locANns & locAN & locBN);% excluding NaN


tsnrBN = tsnrBN(:,locANns & locAN & locBN);
tsnrAN = tsnrAN(:,locANns & locAN & locBN);
tsnrANns = tsnrANns(:,locANns & locAN & locBN);

roi_large = roi_large(locANns & locAN & locBN);

[roi, roiID] = IntersectROIs(roi,roiID,roi_large);  %Is this still necessary?

%% Load BV results for crosscheck that we know how to compute the betas comparing with BV no serial autocorrelation
if(designCheck == 1)
    selrun = 1;
    X  = xff([datafolder,subj,filesep,'Prt_files',filesep,subj,'_run',num2str(selrun),'_Cut_nogap_Collapsed.prt']);
    X1 = design_matrix_BV(X,VTCInfo.NrOfVolumes, VTCInfo.TR);
    % glmBN       = xff([datafolder,filesep,subj,filesep,'GLM_forBVcompare',filesep,subj,'_run',num2str(selrun),'_forBVcompare.glm' ]);
    glmBN       = xff([datafolder,filesep,subj,filesep,'GLM_forBVcompare',filesep,subj,'_run',num2str(selrun),'_forBVcompare_noserialcorr.glm' ]);
    Xselrun     = glmBN.DesignMatrix;
    figure,plot(Xselrun(:,1));hold all;plot(X1(:,1))
    BselrunA     = glmBN.GLMData.BetaMaps;
    tsbncheck    = tsbn; % this is the transformation: we should if standardize with the mean of all runs or the mean run
    Bhat         = inv(Xselrun'*Xselrun)*Xselrun'*tsbncheck(1:230,:);% only for run1 valid
    BselrunROIA  = permute( double(BselrunA) , [4 1 2 3] );  BselrunROIA  = BselrunROIA(:,roi);
    BselrunROIA = BselrunROIA(:,locANns & locAN & locBN);
    figure('Color',[1,1,1]); plot(BselrunROIA(:),Bhat(:),'.');ylim([min(BselrunROIA(:)) max(BselrunROIA(:))])
    xlabel('BV');ylabel('Matlab'); %hold all;plot(BselrunROIA(:),BselrunROIA(:),'r.');
end
%% GLM: general linear model
H           = inv(Xbdiag'*Xbdiag)*Xbdiag';         % no column in X should be only zeros
bBn = H*tsbn;  bAn = H*tsan;   bAnNS = H*tsanNS;   % betas all methods: before after and no noise
ccInv = inv(Xbdiag'*Xbdiag);
%% create maps for BV (large ROI) and saves them all
% F-Maps
connames             = {'Comp_H' , 'Comp_L', 'Oddball_L','Oddball_H', 'Unexp_H', 'Unexp_L' };% names of the contrasts
[Fbn,Fan,FanNS] = SaveFMaps(connames,colnames,Xbdiag,bBn,bAn,bAnNS,tsbn,tsan,tsanNS,ccInv,folderOut,roi_large,VTCInfo,subj);

% single Beta Maps
SaveBMaps(SaveBetaMapsMode,trial_x_run,connames,colnames,bBn,bAn,bAnNS,folderOut,roi_large,VTCInfo,subj);

thF =tinv(1-0.05,size(Xbdiag,1)-size(Xbdiag,2));  %% 0.05 one sided
% tonotopies
connames             = {'Comp_H' , 'Comp_L'; 'Oddball_H','Oddball_L'; 'Unexp_H', 'Unexp_L'}; % names of the contrasts
SaveTono(thF,connames,colnames,bBn,bAn,bAnNS,Fbn,Fan,FanNS,folderOut,roi_large,VTCInfo,subj);

%%
resN   = tsbn - tsan   ;                  % nordic residuals:  Before - After
resNns = tsbn - tsanNS ;                  % nordic residuals:  Before - After no noise
bresN   = H*resN;                         % betas of the residuals
bresNns = H*resNns;                       % betas of the residuals

SaveTonoRes(thF,connames,colnames,bresN,bresNns,Fbn,folderOut,roi_large,VTCInfo,subj);

%% cut to small ROI
tsbn = tsbn(:,roi);
tsan = tsan(:,roi);
tsanNS = tsanNS(:,roi);
bBn = bBn(:,roi);
bAn = bAn(:,roi);
bAnNS = bAnNS(:,roi);
tsnrBN = tsnrBN(:,roi);
tsnrAN = tsnrAN(:,roi);
tsnrANns = tsnrANns(:,roi);
%%
[CorrBetas] = CorrBetasAcrossCond(connames,colnames,bBn,bAn,bAnNS);

groupCorrBetas = 0.*CorrBetas; groupCorrBetas(1,:,:) = 1; groupCorrBetas(2,:,:) = 2; 
groupCorrBetas2 = 0.*CorrBetas; groupCorrBetas2(:,1,:) = 1; groupCorrBetas2(:,2,:) = 2; groupCorrBetas2(:,3,:) = 3; groupCorrBetas2(:,4,:) = 4; 
groupCorrBetas2(:,5,:) = 5; groupCorrBetas2(:,6,:) = 6; 


figure('Color',[1,1,1],'Position',[281 19 1160 820]),

mb = boxchart(groupCorrBetas2(:),CorrBetas(:), 'GroupByColor', groupCorrBetas(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');
ylabel('Spatial Correlation To Original Data','interpreter','Latex')
xlabel('Condition','interpreter','Latex')
legend({'Nordic', 'Nordic No Noise' })


% figure('Color',[1,1,1],'Position',[281 19 1160 820]),
% plot(CorrBetas(1,:),'s','MarkerSize',10,'MarkerFaceColor',[0.3010 0.7450 0.9330], 'MarkerEdgeColor', [0.3010 0.7450 0.9330]); hold on, plot(CorrBetas(2,:),'o','MarkerSize',10,'MarkerFaceColor', [0.4660, 0.6740, 0.1880], 'MarkerEdgeColor', [0.4660 0.6740 0.1880]);
% legend('Nordic-original', 'Nordic NoNoise-original')
% xlim([0.5 6.5])
% ylim([0.2 1])
% xlabel('Condition','interpreter','Latex') 
% ylabel('Correlation','interpreter','Latex')
% title('Spatial Correlation Before-After Nordic','interpreter','Latex')
% 
% s = gcf;
% saveas(s,[folderOut, [subj,'_SpatialCorrelation_BetasAcrossCond.fig']]);
% save([folderOut,[subj, '_SpatialCorrelation_BetasAcrossCond.mat']],'CorrBetas');
% saveas(s, [folderOut, [subj,'_SpatialCorrelation_BetasAcrossCond.svg']]);

%%
resN   = tsbn - tsan   ;                  % nordic residuals:  Before - After
resNns = tsbn - tsanNS ;                  % nordic residuals:  Before - After no noise
tvar  = var(tsbn');                       % checking the variance in the time domain
vvar  = var(tsbn);                        % variance across voxels
bresN   = H*resN;                         % betas of the residuals
bresNns = H*resNns;                       % betas of the residuals
%% here we compute the information regarding the design matrix before nordic and after and in the residuals - leaking

TSSbn = sum( ( tsbn - mean(tsbn) ).^2  );  TSSan = sum( ( tsan - mean(tsan) ).^2  );  TSSanNS = sum( ( tsanNS - mean(tsanNS) ).^2  );
[TSSbn,  SSbbn,  SSebn]   = partitionSSols(tsbn,  Xbdiag,bBn);   Rbn   = SSbbn./TSSbn;
[TSSan,  SSban,  SSean]   = partitionSSols(tsan,  Xbdiag,bAn);   Ran   = SSban./TSSan;
[TSSanNS,SSbanNS,SSeanNS] = partitionSSols(tsanNS,Xbdiag,bAnNS); RanNS = SSbanNS./TSSanNS;

[TSSres,    SSbres,    SSeres]         = partitionSSols(resN ,  Xbdiag,bresN);       Rbres   = SSbres./TSSres;
[TSSresNS,  SSbresNS,  SSeresNS]       = partitionSSols(resNns ,Xbdiag,bresNns);     RbresNS = SSbresNS./TSSresNS;

figure('Color',[1,1,1],'Position',[281 19 1160 820]),
subplot(221);      mh1 = histogram(Rbn,100,'FaceColor',[1 0.5 0.5] ); hold all;histogram(Ran,100,'FaceColor',[0.5 0.5 1] );
xlabel('design matrix SS ratio','interpreter','Latex')
%legend('Original','Nordic');xlim([0 0.6])
legend('Original $\frac{SS_{\mathbf{X}\beta}}{SS_{Y_{O}}}$','Nordic $\frac{SS_{\mathbf{X}\beta}}{SS_{Y_{NO}}}$', 'interpreter','Latex','FontSize',16);
subplot(222);mh1 = histogram(Ran,100,'FaceColor',[0.5 0.5 1] ); hold all;histogram(Rbres,100,'FaceColor',[0.5 1 0.5] );
%legend('Nordic','residuals'); xlabel('design matrix SS ratio','interpreter','Latex')
legend('Nordic $\frac{SS_{\mathbf{X}\beta}}{SS_{Y_{NO}}}$','residuals $\frac{SS_{\mathbf{X}\beta_{\epsilon}}}{SS_{\epsilon}}$', 'interpreter','Latex','FontSize',16);xlim([0 0.6])
subplot(223);   mh1 = histogram(Rbn,100,'FaceColor',[1 0.5 0.5] ); hold all;histogram(RanNS,100,'FaceColor',[0.5 0.5 1] );
legend('Original','Nordic NoNoise');xlabel('design matrix SS ratio','interpreter','Latex')
xlim([0 0.6])
subplot(224);histogram(RanNS,100,'FaceColor',[0.5 0.5 1] ); hold all;histogram(RbresNS,100,'FaceColor',[0.5 1 0.5] );
legend('Nordic NoNoise','residuals NoNoise');xlabel('design matrix SS ratio','interpreter','Latex')
xlim([0 0.6])

h = gcf;
saveas(h,[folderOut, [subj,'_SumOfS_Model.fig']]);
save([folderOut,[subj, '_SumOfS_Model.mat']],'Rbn','Ran','Rbres','RanNS','RbresNS');
saveas(h, [folderOut, [subj,'_SumOfS_Model.svg']]);

%% Create new histogram of the residuals

Rbres_refSSY = SSbres./TSSbn;
RbresNS_refSSY = SSbresNS./TSSbn;

figure('Color',[1,1,1],'Position',[281 19 1160 820]),
subplot(221);      mh1 = histogram(Rbres_refSSY,100,'FaceColor',[1 0.5 0.5] ); hold all;histogram(RbresNS_refSSY,100,'FaceColor',[0.5 0.5 1] );
xlabel('design matrix SS ratio','interpreter','Latex')
%legend('Original','Nordic');
xlim([0 0.6])
legend('Nordic res $\frac{SS_{\mathbf{X}\beta}}{SS_{Y_{O}}}$','Nordic NoNoise res $\frac{SS_{\mathbf{X}\beta}}{SS_{Y_{O}}}$', 'interpreter','Latex','FontSize',16);
subplot(222);mh1 = histogram(Ran,100,'FaceColor',[0.5 0.5 1] ); hold all;histogram(Rbres,100,'FaceColor',[0.5 1 0.5] );
%legend('Nordic','residuals'); xlabel('design matrix SS ratio','interpreter','Latex')
legend('Nordic $\frac{SS_{\mathbf{X}\beta}}{SS_{Y}}$','residuals $\frac{SS_{\mathbf{X}\beta_{\epsilon}}}{SS_{\epsilon}}$', 'interpreter','Latex','FontSize',16);xlim([0 0.6])
subplot(223);   mh1 = histogram(Rbn,100,'FaceColor',[1 0.5 0.5] ); hold all;histogram(RanNS,100,'FaceColor',[0.5 0.5 1] );
legend('Original','Nordic NoNoise');xlabel('design matrix SS ratio','interpreter','Latex')
xlim([0 0.6])
subplot(224);histogram(RanNS,100,'FaceColor',[0.5 0.5 1] ); hold all;histogram(RbresNS,100,'FaceColor',[0.5 1 0.5] );
legend('Nordic NoNoise','residuals NoNoise');xlabel('design matrix SS ratio','interpreter','Latex')
xlim([0 0.6])

% h = gcf;
% saveas(h,[folderOut, [subj,'_SumOfS_Model.fig']]);
% save([folderOut,[subj, '_SumOfS_Model.mat']],'Rbn','Ran','Rbres','RanNS','RbresNS');
% saveas(h, [folderOut, [subj,'_SumOfS_Model.svg']]);


%% Plot histogram of the residuals --> Look at this with Tin
figure('Color',[1,1,1],'Position',[281 19 1160 820]),

histogram(Rbres,100,'FaceColor',[0.3010 0.7450 0.9330]); hold all;histogram(RbresNS,100,'FaceColor',[0.4660 0.6740 0.1880]);
legend('Nordic residuals','Nordic NoNoise residuals');xlabel('design matrix SS ratio','interpreter','Latex');
ylabel('Number of voxels')
xlim([0 0.6])

p = gcf;
saveas(p,[folderOut, [subj,'_SumOfS_Residuals.fig']]);
save([folderOut,[subj, '_SumOfS_Residuals.mat']],'Rbn','Ran','Rbres','RanNS','RbresNS');
saveas(p, [folderOut, [subj,'_SumOfS_Residuals.svg']]);


%% Put values of histograms in a vmp

Ran_vmp = [Ran*10]';
RanNS_vmp = [RanNS*10]';

create_vmp(folderOut,[subj,'_Histogram.vmp'],VTCInfo,[Ran_vmp, RanNS_vmp],roi_large,{'Ran', 'RanNS'})



%% Betas and t consistency
connames             = {'Comp_H' , 'Comp_L', 'Oddball_L','Oddball_H', 'Unexp_H', 'Unexp_L' };% names of the contrasts
connamesNoUnderscore = {'Comp H' , 'Comp L', 'Oddball L','Oddball H', 'Unexp H', 'Unexp L' };% names for the plot
ROIname    = {'HG','PP','PT', 'aSTG','pSTG'};
splithalf_exp     = zeros(size(colnames));% this vector defines which trials are in the first half of the experiment
splithalf_exp( 1 : length(splithalf_exp)/2)=1;
nrep     = 50;
p0       = 0; % percentile for voxel selection, 0-100 if no selection 0
[tbyCond ,tbyCondSH1,tbyCondSH2,BetasCond ] = betas_t_consistency(colnames,connames,connamesNoUnderscore,bBn,bAn,bAnNS,splithalf_exp,roiID, trial_x_run, nrep,ROIname ,p0,folderOut, subj);
%% Doing tonotopic maps here
tonoconH = {'Comp_H','Oddball_H' ,'Unexp_H'};
tonoconL = {'Comp_L','Oddball_L' ,'Unexp_L'};
ROIname    = {'HG','PP','PT','aSTG','pSTG',};
%% within ROI standardization
p0 = 80;
[tonoBN, tonoAN,rsaBN,rsaAN,rsaANnn,prefBN,prefAN,prefANnn] = get_tono(colnames, tonoconH,tonoconL,tbyCond,bBn,bAn,bAnNS,splithalf_exp,roiID,ROIname,trial_x_run,p0, folderOut);
% save('D:\nordic\scripts\temp.mat','tonoBN','tonoAN','rsaBN');
% whole cortex standardization for tonotopy
%% RSA figure
tonocon = {'Comp','OddB','Unexp'};
figure('Color',[1,1,1],'Position',[281 19 1160 820]),
selroi = 3;
% rsaBN = permute(rsaBN,[3, 1, 2]);rsaAN = permute(rsaAN,[3, 1, 2]);rsaANnn = permute(rsaANnn,[3, 1, 2]);
subplot(2,3,1);imagesc(squeeze(rsaBN(selroi,:,:)));
title([ROIname{selroi} 'selectivity RSA: BN-AN'],'interpreter','Latex');colorbar;caxis([0 1])
h = gca(); set(h,'XTick',1:3);set(h,'YTick',1:3);set(h,'XTickLabels',tonocon);set(h,'YTickLabels',tonocon);

subplot(2,3,2);imagesc(squeeze(rsaAN(selroi,:,:)));title([ROIname{selroi} 'selectivity RSA: BN-ANnn'],'interpreter','Latex');colorbar;caxis([0 1])
h = gca(); set(h,'XTick',1:3);set(h,'YTick',1:3);set(h,'XTickLabels',tonocon);set(h,'YTickLabels',tonocon);

subplot(2,3,3);imagesc(squeeze(rsaANnn(selroi,:,:)));title([ROIname{selroi} 'selectivity RSA: AN-ANnn'],'interpreter','Latex');colorbar;caxis([0 1])
h = gca(); set(h,'XTick',1:3);set(h,'YTick',1:3);set(h,'XTickLabels',tonocon);set(h,'YTickLabels',tonocon);

h = gcf;
saveas(h,[folderOut,[subj,'_TonotopyRSA.fig']]);
saveas(h, [folderOut, [subj,'_TonotopyRSA.svg']]);
save([folderOut,[subj,'_TonotopyRSA.mat']],'tonocon','rsaBN', 'ROIname','selroi', 'rsaAN','rsaANnn');

%% Selectivity
figure('Color',[1,1,1],'Position',[281 19 1160 820]),
subplot(231);plot(prefBN{3,1},prefAN{3,1},'.');  xlabel('Before Nordic','interpreter','Latex');ylabel('After Nordic','interpreter','Latex');
subplot(232);plot(prefBN{3,1},prefANnn{3,1},'.');xlabel('Before Nordic','interpreter','Latex');ylabel('After Nordic NN','interpreter','Latex');
subplot(233);plot(prefAN{3,1},prefANnn{3,1},'.');xlabel('After Nordic','interpreter','Latex'); ylabel('After Nordic NN','interpreter','Latex');

h = gcf;
saveas(h,[folderOut,[subj,'_SelectivityPreference.fig']]);
saveas(h, [folderOut, [subj,'_SelectivityPreference.svg']]);
save([folderOut,[subj, '_SelectivityPreference.mat']],'prefBN','prefAN', 'prefANnn');
%% whole brain nordic residuals:  Before - After
p0 = 80;
[tonoBN, tonoAN] = get_tono_wholeBstd(colnames, tonoconH,tonoconL,tbyCond,tbyCondSH1,tbyCondSH2,bBn,bAn,splithalf_exp,roiID,ROIname,p0,trial_x_run);
%% writtig files to BV
% index = repmat(roi,size(tonoBN,1),1);
% mapname = {'Comp','OddBall','Unexpected'};
% create_vmp(['D:\nordic\vmpchecks\tonocheck.vmp'],VTCInfo,10*tonoBN,index,mapname );
% index = repmat(roi,size(tonoBN(1,:),1),1);
% create_vmp(['D:\nordic\vmpchecks\tmapcheck.vmp'],VTCInfo,mtalltrials(1,:),index);
% index = repmat(roi,size(tonoBN(1,:),1),1);% above t
% create_vmp(['D:\nordic\vmpchecks\tmapcheckselection.vmp'],VTCInfo,10*tonoBN(1,mtalltrials(1,:) > 1.65),index(mtalltrials(1,:) > 1.65));
%% beta histograms and t-histograms
figure('Color',[1,1,1]),
for it = 1:6
    subplot(2,6,it);
    histogram(squeeze(BetasCond(it,1,:)),500); hold all;histogram(squeeze(BetasCond(1,2,:)),500);
    xlim([-10 10]);
    title(['$\beta$-cond ', connamesNoUnderscore{it}],'interpreter','Latex')
    if(it==1) legend('before','after'); end

    subplot(2,6,6+it)
    histogram(squeeze(tbyCond(it,1,:)),100); hold all;histogram(squeeze(tbyCond(1,2,:)),100);
    xlim([-5 10]);
    title(['t-cond ', connamesNoUnderscore{it} ],'interpreter','Latex');
end

h = gcf;
saveas(h,[folderOut,[subj, '_Histograms_BetaANDt.fig']]);
saveas(h, [folderOut, [subj,'_Histograms_BetaANDt.svg']]);
save([folderOut, [subj, 'Histograms_BetaANDt.mat']],'BetasCond','connamesNoUnderscore', 'tbyCond');

%% scatter plots and bdiff vs correlation
figure('Color',[1,1,1]),
for it = 1%:6
    %     subplot(3,6,it);
    subplot(2,2,it);
    plot(squeeze(BetasCond(it,1,:)) , squeeze(BetasCond(it,2,:)),'.'); hold all; xl = xlim();
    plot(xl,xl,'r'); 
    xlabel('Before Nordic','interpreter','Latex'); ylabel('After Nordic','interpreter','Latex')
    hold all;title(['$\beta$-cond ', connamesNoUnderscore{it}],'interpreter','Latex'); xlim(xl)
    xlim([-50 50]);ylim([-50 50]); axis square
    %     subplot(3,6,6+it)
    subplot(2,2,2)
    con = contains(colnames, connames{it});%get the exact columns where each condition occurs
    cResDesign = corr(sum(Xbdiag(:,con),2),resN);%correlations of residuals with design of the particular exp conditions
    plot(squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:)),cResDesign,'.'); hold all; xl = xlim();
    xlabel('Before - after','interpreter','Latex');
    ylabel('corr Nordic res(bef-after) x design','interpreter','Latex');
    title('corr res vs coldesign','interpreter','Latex');xlim([-10 10]) ;axis square
    %     subplot(3,6,12+it)
    subplot(2,2,3)
    histogram(squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:)) , 500); hold all
    mdiff = mean(squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:)));
    plot([mdiff mdiff],ylim(),'LineWidth',2);
    xlabel('Before - after','interpreter','Latex'); xlim([-10 10]) ;axis square
    ylabel('Num voxel','interpreter','Latex');      xlim([-10 10]) ;
    title('hist of $\beta$ diff','interpreter','Latex');

    subplot(2,2,4)
    plot(squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:)) , mean(tsnrBN ) ,'.'); hold all
    [my binsx ] = window_mean(  squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:)),  mean(tsnrBN), linspace(5,80,20) );
    plot(my,binsx,'r-','LineWidth',2)
    xlabel('Before - after','interpreter','Latex');%$ xlim([-10 10]) ;
    axis square
    ylabel('tSNR','interpreter','Latex');      xlim([-10 10]) ;

end

h = gcf;
saveas(h,[folderOut,[subj, '_ScatterplotsBetas_corrwDesign_tSNR.fig']]);
saveas(h, [folderOut, [subj,'_ScatterplotsBetas_corrwDesign_tSNR.svg']]);
save([folderOut, [subj,'_ScatterplotsBetas_corrwDesign_tSNR.mat']],'BetasCond','connamesNoUnderscore', 'con', 'colnames', 'connames', 'Xbdiag', 'cResDesign', 'mdiff', 'tsnrBN', 'my', 'binsx');