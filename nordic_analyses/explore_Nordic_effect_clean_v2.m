clc;clear all;close all
%% set path
datafolder = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\';

SaveBetaMapsMode = 1;    % 1 saves average betas per run per condition - 2 saves single trials (separated per condition in different maps)

subj       = 'S2_MN_NG_PC';     % Change participant number accordingly
nruns      = 8;                 % Change number of runs accordingly --> S5 and S6 have 6 runs, S10 has 7 runs and the rest has 8

% outputfolder
folderOut = [datafolder,subj,filesep,'ResultsScriptVolume',filesep];
if ~exist(folderOut)
    mkdir(folderOut);
end

xff(0,'transiosize','vtc',1048576);% read vtc fast
designCheck = 0;


%% Load ROIs from these names: function extract_roi_indices
roiAll      = {'TemporalLobe.msk'};  %Very large mask of temporal lobe including WM and CSF
roinames    = {'HG_GM.msk','PP_GM.msk','PT_GM.msk','aSTG_GM.msk','pSTG_GM.msk'};  %ROI's confined to GM coming from segmentation
%        roinames    = {'HG_GM.msk'};
[roi_large, ~]   = extract_roi_indices(roiAll,datafolder,subj);
[roi_dup, roiID] = extract_roi_indices(roinames,datafolder,subj);

%% Remove duplicates from roi's
ii = roi_dup == roi_dup';

[roi, IA, IC] = unique(roi_dup, 'first');
roiID = roiID(IA);
%%
p0exclude = 0.5; % in percentage
%% Loads time series data: before NORDIC (bN), after NORDIC (aN), after NORDIC no noise (AnNS)
tsbn     = []; tsan = []; tsanNS = [];
Xbdiag   = [];  colnames = {};
runIdx = [];
locBN    = ~isnan(roi_large);locAN    = ~isnan(roi_large);  locANns = ~isnan(roi_large);
trial_x_run = [];
for itr    = 1:nruns
    disp(['Reading run ' , num2str(itr)]);
    % Load VTC before NORDIC of each run: function read_vtc which also does  % signal change standardization
    [tsbn_t ,VTCInfo,tsnrBN(itr,:)] = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc' ],roi_large);
    %[tsbn_t ,VTCInfo,tsnrBN(itr,:)] = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.vtc' ],roi_large);
    %[tsbn_t ,VTCInfo,tsnrBN(itr,:)] = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c.vtc' ],roi_large);
    %[tsbn_t ,VTCInfo,tsnrBN(itr,:)] = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_Cut_SCSTBL_3DMCTS.vtc' ],roi_large);

    voxExclude           = tsnrBN(itr,:) < prctile(tsnrBN(itr,:),p0exclude );% for every run localize the voxels in the lowest p0exclude percent in tsnrBN
    tsbn_t(:,voxExclude) = NaN; % and make this voxels NaN since they will be systematically excluded in all the analyses
    locBN                           = locBN & sum(isnan(tsbn_t ))==0; % locate NaN in the data
    % load after NORDIC
    [tsan_t,~,  tsnrAN(itr,:)]      = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc' ],roi_large);
    %[tsan_t,~,  tsnrAN(itr,:)]      = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.vtc' ],roi_large);
    %[tsan_t,~,  tsnrAN(itr,:)]      = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c.vtc' ],roi_large);
    %[tsan_t,~,  tsnrAN(itr,:)] = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NORDIC_Cut_SCSTBL_3DMCTS.vtc' ],roi_large);
    locAN             = locAN & sum(isnan(tsan_t ))==0;% locate NaN in the data
    % load after NORDIC no noise
    [tsanNS_t,~,tsnrANns(itr,:)]   = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc' ],roi_large);
    %[tsanNS_t,~,tsnrANns(itr,:)]   = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.vtc' ],roi_large);   
    %[tsanNS_t,~,tsnrANns(itr,:)]   = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c.vtc' ],roi_large);
    %[tsanNS_t,~,tsnrANns(itr,:)]   = read_vtc([datafolder,subj,filesep,'VTC',filesep,subj,'_run',num2str(itr),'_NoN_Cut_SCSTBL_3DMCTS.vtc' ],roi_large);
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
tsbn   = tsbn(:,  locANns & locAN & locBN);  % time series before NORDIC, excluding NaN
tsan   = tsan(:,  locANns & locAN & locBN);  % excluding NaN
tsanNS = tsanNS(:,locANns & locAN & locBN);% excluding NaN


tsnrBN = tsnrBN(:,    locANns & locAN & locBN); % time series per run.
tsnrAN = tsnrAN(:,    locANns & locAN & locBN);
tsnrANns = tsnrANns(:,locANns & locAN & locBN);

roi_large = roi_large(locANns & locAN & locBN);

[roi, roiID] = IntersectROIs(roi,roiID,roi_large);  

%% GLM: general linear model
H     = inv(Xbdiag'*Xbdiag)*Xbdiag';         % no column in X should be only zeros
bBn   = H*tsbn;  bAn = H*tsan;   bAnNS = H*tsanNS;   % betas all methods: before after and no noise
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

%% Calculate residuals for NORDIC
% resN   = tsbn - tsan   ;                  % NORDIC residuals:  time series before NORDIC - time series After NORDIC
% resNns = tsbn - tsanNS ;                  % NORDIC residuals:  Before - After no noise
mbn    = mean(tsbn); man = mean(tsan);      % take the mean from both time series --> in PSC so mean is 100
tsbnC       = tsbn - mbn; tsanC       = tsan - man;  % Centered time series, mean = 0.
alphaN1     = sum(tsbnC.*tsanC)./sum(tsanC.^2);     % finding the scaling factor alpha
% checking we are doing OLS fine against GLMfit for a random voxel
% voxsel      = 1023;
% yhatL       = alphaN1(voxsel)*(tsan(:,voxsel) - man(voxsel)) + mbn(voxsel);% predictions of my GLM model
% alphaCheck = glmfit(tsan(:,voxsel) ,tsbn(:,voxsel) ); % OLS in matlab
% yhatM      = alphaCheck(1) + alphaCheck(2).*tsan(:,voxsel); % predictions OLS matlab
% figure('Color',[1,1,1]);plot(tsan(:,voxsel),tsbn(:,voxsel),'.'); hold all% cross check
% plot(tsan(:,voxsel) , yhatL,'.');plot(tsan(:,voxsel) , yhatM,'.')
resN = tsbn - alphaN1.*(tsan - man) - mbn;          % New definition of residuals, mean of residuals = 0
yR   =  alphaN1.*(tsan - man) + mbn + resN;         
ssYR = sum( (alphaN1.*(tsan - man)       + mbn    + resN).^2 ) - size(resN,1)*mbn.^2;
ssYR = sum( (alphaN1.^2.*(tsan - man).^2 + mbn.^2 + resN.^2)) + 2*sum(alphaN1.*(tsan - man).*resN) + 2*sum(mbn.*resN) + 2*sum(mbn.*alphaN1.*(tsan - man) )               - size(resN,1)*mbn.^2;
ssYR = sum( (alphaN1.^2.*(tsan - man).^2 + mbn.^2 + resN.^2)) + 0                                  + 0                + 0                                                 - size(resN,1)*mbn.^2;
ssYR = alphaN1.^2.*sum( (tsan - man).^2) + sum(resN.^2);

%% Calculate residuals for NORDIC NN

manNS = mean(tsanNS);      % take the mean from time series --> in PSC so mean is 100
tsanNSC       = tsanNS - manNS;  % Centered time series, mean = 0.
alphaNS1     = sum(tsbnC.*tsanNSC)./sum(tsanNSC.^2);     % finding the scaling factor alpha
% checking we are doing OLS fine against GLMfit for a random voxel
%voxsel      = 1023;

resNns = tsbn - alphaNS1.*(tsanNS - manNS) - mbn;          % New definition of residuals, mean of residuals = 0
yRNS   =  alphaNS1.*(tsanNS - manNS) + mbn + resNns;         % R? is BN or residuals?
ssYRNS = sum( (alphaNS1.*(tsanNS - manNS)       + mbn    + resNns).^2 ) - size(resNns,1)*mbn.^2;
ssYRNS = sum( (alphaNS1.^2.*(tsanNS - manNS).^2 + mbn.^2 + resNns.^2)) + 2*sum(alphaNS1.*(tsanNS - manNS).*resNns) + 2*sum(mbn.*resNns) + 2*sum(mbn.*alphaNS1.*(tsanNS - manNS) )               - size(resNns,1)*mbn.^2;
ssYRNS = sum( (alphaNS1.^2.*(tsanNS - manNS).^2 + mbn.^2 + resNns.^2)) + 0                                  + 0                + 0                                                 - size(resNns,1)*mbn.^2;
ssYRNS = alphaNS1.^2.*sum( (tsanNS - manNS).^2) + sum(resNns.^2);


bresN   = H*resN;                         % betas of the residuals after NORDIC
bresNns = H*resNns;                       % betas of the residuals after NORDIC NN

SaveTonoRes(thF,connames,colnames,bresN,bresNns,Fbn,folderOut,roi_large,VTCInfo,subj);

%% Spatial Correlation of betas across the number of runs per condition (WHOLE BRAIN NO Threshold) --> figures 1, 2 and 3
p0 = 1;

[CorrBetas, tonoBn, tonoAn, tonoAnNS] = CorrBetasAcrossCond(connames,colnames,bBn,bAn,bAnNS,p0,folderOut,roi_large,VTCInfo,subj);


%% Create figure of Spatial Correlations
groupCorrBetas = 0.*CorrBetas;  groupCorrBetas(1,:,:)  = 1; groupCorrBetas(2,:,:) = 2;
groupCorrBetas2 = 0.*CorrBetas; groupCorrBetas2(:,1,:) = 1; groupCorrBetas2(:,2,:) = 2; groupCorrBetas2(:,3,:) = 3; groupCorrBetas2(:,4,:) = 4;
groupCorrBetas2(:,5,:) = 5; groupCorrBetas2(:,6,:) = 6;

figure('Color',[1,1,1],'Position',[281 19 1160 820]),

mb = boxchart(groupCorrBetas2(:),CorrBetas(:), 'GroupByColor', groupCorrBetas(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');
mb(1).BoxFaceColor = [144 103 167]./256;
mb(2).BoxFaceColor = [0.4660 0.6740 0.1880];
ylabel('Spatial Correlation', 'FontName', 'Arial', 'FontSize', 12)
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 12)
xticks([1, 2, 3, 4, 5, 6])
xticklabels({'PredH', 'MispredH', 'UnpredH', 'PredL', 'MispredL', 'UnpredL'})
legend({'NORdef vs Original', 'NORnn vs Original'}, 'FontName', 'Arial', 'FontSize', 10)
ylim([0.3 1])

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';

g = gcf;
% do t-ttest here first all cond together
% [H,PcollapsedCond,CI,STATS] = ttest2(vec( CorrBetas(1,:,: ) ) , vec( CorrBetas(2,:,: ) ) ); % comparing nordic_standard vs nordicNN_vs_Starndard
%           STATcollapsedCond = STATS.tstat;
% % do t-ttest here first all cond together
% [H,P_x_Cond,CI,STATS] = ttest2( squeeze(CorrBetas(1,:,: ) )' , squeeze( CorrBetas(2,:,: ) )' ); % comparing nordic_standard vs nordicNN_vs_Starndard
%           STAT_x_Cond = STATS.tstat;

saveas(g,[folderOut, [subj,'_SpatialCorrelation_BetasAcrossCond.fig']]);
save([folderOut,[subj, '_SpatialCorrelation_BetasAcrossCond.mat']],'CorrBetas');
saveas(g, [folderOut, [subj,'_SpatialCorrelation_BetasAcrossCond.svg']]);



%% here we compute the information regarding the design matrix before NORDIC and after and in the residuals - leaking
TSSbn = sum( ( tsbn - mean(tsbn) ).^2  );  TSSan = sum( ( tsan - mean(tsan) ).^2  );  TSSanNS = sum( ( tsanNS - mean(tsanNS) ).^2  );
[TSSbn,  SSbbn,  SSebn]   = partitionSSols(tsbn,  Xbdiag,bBn);   %Total sum of squares, sum of squares betas and sum of squares error
Rbn   = SSbbn./TSSbn;    % Ratio of sum of squares --> SS of Bn ./ Totals SS of Bn

[TSSan,  SSban,  SSean]   = partitionSSols(tsan,  Xbdiag,bAn);   
Ran   = SSban./TSSan;

[TSSanNS,SSbanNS,SSeanNS] = partitionSSols(tsanNS,Xbdiag,bAnNS); 
RanNS = SSbanNS./TSSanNS;


[TSSres,    SSbres,    SSeres]         = partitionSSols(resN ,  Xbdiag,bresN);    % Same as above but for residuals   
Rbres   = SSbres./TSSres;   % Normalized to TSS in residuals

[TSSresNS,  SSbresNS,  SSeresNS]       = partitionSSols(resNns ,Xbdiag,bresNns);     
RbresNS = SSbresNS./TSSresNS;

Rbres_refSSY = SSbres./TSSbn;    % Normalized SS in residuals to the TSS of BEFORE NORDIC
RbresNS_refSSY = SSbresNS./TSSbn;
%% Understanding variance partitioning 
% voxsel = 12003;
% covYnordicResNordic = 2*( sum(resN.*tsan) - size(tsan,1)*mean(resN).*mean(tsan) ); % This term is currenlty 0.
% SScheck   = (alphaN1.^2).*TSSan  + TSSres + covYnordicResNordic;
% SScheck1  = (alphaN1.^2).*TSSan  + TSSres;
% %  sspart   = [TSSbn(voxsel), (alphaN1(voxsel).^2).*TSSan(voxsel), TSSres(voxsel) SScheck1(voxsel)];
% sspart   = [[0, (alphaN1(voxsel).^2).*SSban(voxsel),SSbres(voxsel), 0 0]; [0, (alphaN1(voxsel).^2).*SSean(voxsel) ,SSeres(voxsel), 0 0] ; [TSSbn(voxsel), 0,0, SScheck1(voxsel) , SScheck(voxsel)] ]' ; 
% barnames = {'TSSbn','TSSan','TSSresN','TSScheck'};
% figure('Color',[1,1,1]);
% bar(sspart(1:4,:),'stacked');
% set(gca,'xticklabel',barnames);
% h = gca();set(h,'FontSize',16)
% title(['Variance partitioning voxel ' , num2str(voxsel)]);
% legend('Design Related','Design Unrelated')
% 
% %% Understanding variance partioning NORDIC NN
% covYnordicResNordicNS = 2*( sum(resNns.*tsanNS) - size(tsanNS,1)*mean(resNns).*mean(tsanNS) ); % This term is currenlty 0.
% SScheck   = (alphaNS1.^2).*TSSanNS  + TSSresNS + covYnordicResNordicNS;
% SScheck1  = (alphaNS1.^2).*TSSanNS  + TSSresNS;
% %  sspart   = [TSSbn(voxsel), (alphaN1(voxsel).^2).*TSSan(voxsel), TSSres(voxsel) SScheck1(voxsel)];
% sspart   = [[0, (alphaNS1(voxsel).^2).*SSbanNS(voxsel),SSbresNS(voxsel), 0 0]; [0, (alphaNS1(voxsel).^2).*SSeanNS(voxsel) ,SSeresNS(voxsel), 0 0] ; [TSSbn(voxsel), 0,0, SScheck1(voxsel) , SScheck(voxsel)] ]' ; 
% barnames = {'TSSbn','TSSanNS','TSSresNns','TSScheck'};
% figure('Color',[1,1,1]);
% bar(sspart(1:4,:),'stacked');
% set(gca,'xticklabel',barnames);
% h = gca();set(h,'FontSize',16)
% title(['Variance partitioning voxel ' , num2str(voxsel)]);
% legend('Design Related','Design Unrelated')
%% Create vmp's of SS of what is in the residuals in relation to variance in Standard data or residuals

Ran_vmp = [Rbres_refSSY*10]';
RanNS_vmp = [RbresNS_refSSY*10]';
create_vmp(folderOut,[subj,'_SSresiduals.vmp'],VTCInfo,[Ran_vmp, RanNS_vmp],roi_large,{'SS residuals NORDIC', 'SS residuals NORDIC NN'})

%% cut to small ROI s
tsbn     = tsbn(:,roi);
tsan     = tsan(:,roi);
tsanNS   = tsanNS(:,roi);
bBn      = bBn(:,roi);
bAn      = bAn(:,roi);
bAnNS    = bAnNS(:,roi);
tsnrBN   = tsnrBN(:,roi);
tsnrAN   = tsnrAN(:,roi);
tsnrANns = tsnrANns(:,roi);

%% Betas and t consistency --> per ROI
connames             = {'Comp_H' , 'Comp_L', 'Oddball_L','Oddball_H', 'Unexp_H', 'Unexp_L' };% names of the contrasts
connamesNoUnderscore = {'Comp H' , 'Comp L', 'Oddball L','Oddball H', 'Unexp H', 'Unexp L' };% names for the plot
ROIname    = {'HG','PP','PT', 'aSTG','pSTG'};
splithalf_exp     = zeros(size(colnames));% this vector defines which trials are in the first half of the experiment
splithalf_exp( 1 : length(splithalf_exp)/2)=1;
nrep     = 50;
p0       = 0; % percentile for voxel selection, 0-100 if no selection 0
[tbyCond ,tbyCondSH1,tbyCondSH2,BetasCond] = betas_t_consistency_new(colnames,connames,connamesNoUnderscore,bBn,bAn,bAnNS,splithalf_exp,roiID, trial_x_run, nrep,ROIname ,p0,folderOut, subj);
                                                

%%
figure('Color',[1,1,1],'Position',[281 19 1160 820]),       % NORDIC vs Standard --> variance explained by design
subplot(121);      mh1 = histogram(Rbn,100); hold all;
histogram(Ran,100);
histogram(RanNS,100);
xlabel('$\frac{Design\ SS}{Total\ SS}$', 'interpreter','Latex','FontSize',20)
ylabel('Number of voxels', 'FontName', 'Arial', 'FontSize', 12)
legend({'Original','NORdef','NORnn'}, 'FontName', 'Arial', 'FontSize', 10);
xlim([0  prctile([Rbn Ran RanNS],99.99)])

% b = gca;
% b.XAxis.FontSize = 12;
% b.XAxis.FontName = 'Arial';
% b.YAxis.FontSize = 12;
% b.YAxis.FontName = 'Arial';

% How much variance explained by design is in the residual time series in
% relation to total SS of the data
subplot(122);      
mh1 = histogram(Rbres_refSSY,linspace(0,1,100),'FaceColor',[144 103 167]./256 ); 
hold all;histogram(RbresNS_refSSY,linspace(0,1,100),'FaceColor',[0.4660 0.6740 0.1880] );
xlabel('$\frac{Design\ SS\ Residuals}{Total\ SS\ Original}$', 'interpreter','Latex','FontSize',20)
ylabel('Number of voxels', 'FontName', 'Arial', 'FontSize', 12)
legend({'NORdef vs Original', 'NORnn vs Original'}, 'FontName', 'Arial', 'FontSize', 10);
xlim([0 0.6])

% b = gca;
% b.XAxis.FontSize = 12;
% b.XAxis.FontName = 'Arial';
% b.YAxis.FontSize = 12;
% b.YAxis.FontName = 'Arial';

% subplot(233);
% it = 1;
% temp          = squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:));
% [Z bins]      = hist3([mean(tsnrBN)',temp],[100,100]); 
% [xgrid ygrid] = meshgrid(bins{2}, bins{1} );
% myp = plot(squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:)) , mean(tsnrBN ) ,'.'); hold all
% myp.Color = [0.301 0.75 0.933];
% [~, myc] = contour(xgrid,ygrid,Z,3, 'k-','LineWidth',2);xlim([-10 10]);axis square
% myc.LineColor = [0.3 0.3 0.3];
% [my binsx ] = window_mean(  squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:)),  mean(tsnrBN), linspace(5,80,10) );
% plot(my,binsx,'r-','LineWidth',2)
% xlabel('Betas Standard - NORDIC', 'FontName', 'Arial', 'FontSize', 12);
% axis square;ylabel('mean/std', 'FontName', 'Arial', 'FontSize', 12);      
% xlim([-10 10]);
% 
% b = gca;
% b.XAxis.FontSize = 12;
% b.XAxis.FontName = 'Arial';
% b.YAxis.FontSize = 12;
% b.YAxis.FontName = 'Arial';

h = gcf;
saveas(h,[folderOut, [subj,'_SumOfS_Model.fig']]);
save([folderOut,[subj, '_SumOfS_Model.mat']],'Rbn','Ran','Rbres','RanNS','RbresNS', 'Rbres_refSSY', 'RbresNS_refSSY', 'tsnrBN', 'BetasCond');
saveas(h, [folderOut, [subj,'_SumOfS_Model.svg']]);


figure()
subplot(121);
it = 1;
temp          = squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:));
[Z bins]      = hist3([mean(tsnrBN)',temp],[100,100]); 
[xgrid ygrid] = meshgrid(bins{2}, bins{1} );
myp = plot(squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:)) , mean(tsnrBN ) ,'.'); hold all
myp.Color = [0.301 0.75 0.933];
% [~, myc] = contour(xgrid,ygrid,Z,3, 'k-','LineWidth',2);xlim([-10 10]);axis square
% myc.LineColor = [0.3 0.3 0.3];
[my binsx ] = window_mean(  squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,2,:)),  mean(tsnrBN), linspace(5,80,10) );
plot(my,binsx,'r-','LineWidth',2)
xlabel('Betas Ori - NORdef', 'FontName', 'Arial', 'FontSize', 12);
axis square;
ylabel('$\frac{Mean}{Std}$', 'interpreter','Latex','fontsize', 20)
xlim([-10 10]);

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';


subplot(122);
it = 1;
temp          = squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,3,:));
[Z bins]      = hist3([mean(tsnrBN)',temp],[100,100]); 
[xgrid ygrid] = meshgrid(bins{2}, bins{1} );
myp = plot(squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,3,:)) , mean(tsnrBN ) ,'.'); hold all
myp.Color = [0.301 0.75 0.933];
% [~, myc] = contour(xgrid,ygrid,Z,3, 'k-','LineWidth',2);xlim([-10 10]);axis square
% myc.LineColor = [0.3 0.3 0.3];
[my binsx ] = window_mean(  squeeze(BetasCond(it,1,:)) - squeeze(BetasCond(it,3,:)),  mean(tsnrBN), linspace(5,80,10) );
plot(my,binsx,'r-','LineWidth',2)
xlabel('Betas Ori - NORnn', 'FontName', 'Arial', 'FontSize', 12);
axis square;
ylabel('$\frac{Mean}{Std}$', 'interpreter','Latex','fontsize', 20)
xlim([-10 10]);

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';

h = gcf;
saveas(h,[folderOut, [subj,'_tSNR.fig']]);
save([folderOut,[subj, '_tSNR.mat']],'tsnrBN', 'BetasCond');
saveas(h, [folderOut, [subj,'_tSNR.svg']]);

