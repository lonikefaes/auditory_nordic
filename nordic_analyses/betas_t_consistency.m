function [tbyCond,tbyCondSH1,tbyCondSH2,BetasCond ] = betas_t_consistency(colnames,connames,connamesNoUnderscore,bBn,bAn,bAnNS,splithalf_exp , roiID, trial_x_run, nrep,roiNames,p0,folderOut,subj)
% Do the reliability analysis of betas and t 
% bBn betas before nordic
% bAn betas after nordic
% connames             = {'Comp_H' , 'Comp_L', 'Oddball_L','Oddball_H', 'Unexp_H', 'Unexp_L' };
% connamesNoUnderscore = {'Comp H' , 'Comp L', 'Oddball L','Oddball H', 'Unexp H', 'Unexp L' };
% splithalf_exp     = zeros(size(colnames));% this vector defines which trials are in the first half of the experiment
% splithalf_exp( 1 : length(splithalf_exp)/2)=1;
colors = [0 0 0.5; 0.3010 0.7450 0.9330; 0.4660, 0.6740, 0.1880];

for itcon = 1:length(connames )
    con = contains(colnames, connames{itcon});
    disp(['doing contrast ', connames{itcon},' involving ', num2str(sum(con)), ' single trials' ]);
    [H,P,CI,STATSbn]      = ttest(bBn(con,:)  ); %tbn(itcon,:) = STATSbn.tstat;
    [H,P,CI,STATSan]      = ttest(bAn(con,:)  ); %tan(itcon,:) = STATSan.tstat;
    [H,P,CI,STATSanNS]    = ttest(bAnNS(con,:)  ); %tan(itcon,:) = STATSan.tstat;
    BetasCond(itcon,:,:)  = [mean(bBn(con,:)) ; mean(bAn(con,:)); mean(bAnNS(con,:)) ];% mean of the single trial betas for each contrast defined by the names in connames
    tbyCond(itcon,:,:)    = [STATSbn.tstat ; STATSan.tstat; STATSanNS.tstat] ; 
    % here we compute the split half correlations of the betas and the t
    shBetas(itcon,:)      =  [ corr( mean(bBn(con & splithalf_exp,:))', mean(bBn(con & ~splithalf_exp,:))'  )          corr( mean(bAn(con & splithalf_exp,:))', mean(bAn(con & ~splithalf_exp,:))'  )    ];
    % split half t maps
     [H,P,CI,STATSbnSH1]  = ttest(bBn(con & splithalf_exp,:)  );     [H,P,CI,STATSbnSH2] = ttest(bBn(con & ~splithalf_exp,:)  );
     [H,P,CI,STATSanSH1]  = ttest(bAn(con & splithalf_exp,:)  );     [H,P,CI,STATSanSH2] = ttest(bAn(con & ~splithalf_exp,:)  );     
     shTmaps(itcon,:)         =  [corr(STATSbnSH1.tstat',STATSbnSH2.tstat'),corr(STATSanSH1.tstat',STATSanSH2.tstat')  ];
     tbyCondSH1(itcon,:,:)    = [STATSbnSH1.tstat ; STATSanSH1.tstat] ; 
     tbyCondSH2(itcon,:,:)    = [STATSbnSH2.tstat ; STATSanSH2.tstat] ; 
     BetasCondSH1(itcon,:)   =  mean( bBn(con & splithalf_exp,:) );
     BetasCondSH2(itcon,:)   =  mean( bBn(con & ~splithalf_exp,:) ); 

     
     
end

mBetasCond = squeeze(mean(BetasCond,3)); % mean across voxels inside ROI of betas for each cond before and after
mtbyCond   = squeeze(mean(tbyCond,3)); %
sdCond     = squeeze(std(BetasCond,[],3));
nvoxROI    = size(bBn,2);
%% 
uroi = unique(roiID);
trial_x_run_dummy = dummyvar(trial_x_run);
nruns             = size(trial_x_run_dummy,2);
conCompL          = contains(colnames, connames{1});
conCompH          = contains(colnames, connames{2});
con               =  conCompL | conCompH;% Complete H Complete L
for itrep = 1:nrep
    runsperm       = randperm(nruns);
    splithalf_exp  = sum(trial_x_run_dummy(:,runsperm(1:nruns/2)),2)';
%     %no split variability
            mbBsh1 = mean(bBn(con &  splithalf_exp,:));
            mbBsh2 = mean(bBn(con & ~splithalf_exp,:));
            mbAsh1 = mean(bAn(con &  splithalf_exp,:));
            mbAsh2 = mean(bAn(con & ~splithalf_exp,:)); 
            mbAnssh1 = mean(bAnNS(con &  splithalf_exp,:));
            mbAnssh2 = mean(bAnNS(con & ~splithalf_exp,:)); 
 
           [H,P,CI,STATSbnSH1]    = ttest(bBn(con & splithalf_exp,:)  );     [H,P,CI,STATSbnSH2] = ttest(bBn(con & ~splithalf_exp,:)  );
           [H,P,CI,STATSanSH1]    = ttest(bAn(con & splithalf_exp,:)  );     [H,P,CI,STATSanSH2] = ttest(bAn(con & ~splithalf_exp,:)  );     
           [H,P,CI,STATSanNSSH1]  = ttest(bAnNS(con & splithalf_exp,:)  );     [H,P,CI,STATSanNSSH2] = ttest(bAnNS(con & ~splithalf_exp,:)  );     


    for itroi = 1:length(uroi)
        locroi = find(roiID  == uroi(itroi));
        if(p0 > 0)  % select those voxels above a percentile thres at each ROI
        selvox = STATSbnSH1.tstat( locroi ) > prctile( STATSbnSH1.tstat( locroi ) , p0);
        locroi = locroi(selvox);
        end; 
     relBetas(itroi,itrep,:) = [ corr(mbBsh1(locroi)' ,mbBsh2(locroi)' )  corr(mbAsh1(locroi)' ,mbAsh2(locroi)' )  corr(mbAnssh1(locroi)' ,mbAnssh2(locroi)' )  ];
     relTs(itroi,itrep,:)    = [ corr(STATSbnSH1.tstat(locroi)' ,STATSbnSH2.tstat(locroi)' )  corr(STATSanSH1.tstat(locroi)' ,STATSanSH2.tstat(locroi)' )  corr(STATSanNSSH1.tstat(locroi)' ,STATSanNSSH2.tstat(locroi)' )];
     roiCount(itroi,itrep,:) = [itroi,itroi,itroi] ;
     methCount(itroi,itrep,:)= [1 2 3];
     BetasCondrepSH(itroi,itrep,:) =  [median( [  mbBsh2(locroi)' ])      median( mbAsh2(locroi)')       median( mbAnssh2(locroi)')  ];
     tCondrepSH(itroi,itrep,:)     =  [median(  STATSbnSH2.tstat(locroi)')           median(STATSanSH2.tstat(locroi)')   median(STATSanNSSH2.tstat(locroi)')   ];
    end
end
%% plot
figure('Color',[1,1,1],'Position',[81 59 2060 820]),
mBetasCompleteLH = squeeze( mean( BetasCond([1,2],:,:) ));% Complete L and H
meth             = 0.*mBetasCompleteLH; meth(1,:) = 1; meth(2,:)=2;meth(3,:)=3;
rois             = 2.*repmat(roiID,size(mBetasCompleteLH,1),1);
%  subplot(221);
subplot(221); mb = boxchart(rois(:),mBetasCompleteLH(:),'GroupByColor', meth(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');
colororder([colors]);
%subplot(221);mb = boxchart(rois(:),mBetasCompleteLH(:),'GroupByColor', meth(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');

%   ylim([-200 400])
  ylim([-2 6])
uroi = unique(roiID);
h   = gca; set(h,'FontSize',12);
set(h,'XTick',unique(rois(:)))
set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
title('$\beta$ change',    'interpreter','Latex')
ylabel('$\beta$  \% change','interpreter','Latex')

 %%
% hold all;
% tt = squeeze(mean(BetasCondrepSH,2));
% plot( 2*squeeze(roiCount(:,1,1)) , tt(:,1),'.','MarkerSize',20)
% for itroi = 1:length(uroi)
%         locroi = find(roiID  == uroi(itroi));
%         mBetas1(itroi) = mean( mean( BetasCondSH1(:,locroi)) );
%         mBetas2(itroi) = mean( mean( BetasCondSH2(:,locroi)) );
% end
% plot( 2*squeeze(roiCount(:,1,1)) , mBetas1,'bo','MarkerSize',10)
% plot( 2*squeeze(roiCount(:,1,1)) , mBetas2,'g.','MarkerSize',20)

% subplot(222)
% mb = boxchart(2.*roiCount(:),BetasCondrepSH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
% ylim([-1 4])

%%
rois = 2.*uroi;
figure('Color',[1,1,1],'Position',[81 59 2060 820]),
subplot(221)
boxchart(2.*roiCount(:),BetasCondrepSH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
colororder([colors]);

% mb  = bar(mBetasCond); 
% legend('before','after'); title('$\beta$','interpreter','Latex');
% ylim([prctile(mBetasCompleteLH(:),10),  prctile(mBetasCompleteLH(:),90)])
ylim([-1.5 7]);
h   = gca; set(h,'FontSize',12); 
set(h,'XTick',unique(rois(:)))
set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
title('$\beta$ change',    'interpreter','Latex')
ylabel('$\beta$  \% change','interpreter','Latex')

subplot(222);
mtCompleteLH = squeeze(mean(tbyCond(1:2,:,:),1));
boxchart(2.*roiCount(:),tCondrepSH(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
colororder([colors]);
title('t-values','interpreter','Latex'); hold all;
h = gca; set(h,'FontSize',12);
% ylim([prctile(tCondrepSH(:),1), prctile(tCondrepSH(:),99)])
ylim([-1 8]);
set(h,'XTick',unique(rois(:)))
set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
title('$t$ change','interpreter','Latex')
ylabel('t-statistic','interpreter','Latex')

% errorbar(mb(1,1).XEndPoints,mtbyCond(:,1),sdCond(:,1),'LineWidth',2)
subplot(223);
boxchart(2.*roiCount(:),relBetas(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
colororder([colors]);

title('$\beta$ reliability','interpreter','Latex'); hold all;
h = gca; set(h,'FontSize',12);
ylim([-0.2 1])
set(h,'XTick',unique(rois(:)))
set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
ylabel('split half corr','interpreter','Latex')
% bar( shBetas);title('$\beta$ split half corr (across vox)','interpreter','Latex'); hold all;
% h = gca; set(h,'FontSize',12);set(h,'XTickLabel',connamesNoUnderscore,'TickLabelInterpreter','Latex')

subplot(224);
boxchart(2.*roiCount(:),relTs(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
colororder([colors]);

title('t-value reliability','interpreter','Latex'); hold all;
h = gca; set(h,'FontSize',12);
ylim([-0.2 1])
set(h,'XTick',unique(rois(:)))
set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
ylabel('split half corr','interpreter','Latex')

h = gcf;
saveas(h,[folderOut, [subj,'_BetasANDt_changeANDReliability.fig']]);
saveas(h, [folderOut, [subj,'_BetasANDt_changeANDReliability.svg']]);
save([folderOut,[subj, '_BetasANDt_changeANDReliability.mat']],'uroi','roiCount','BetasCondrepSH','methCount','relBetas','relTs','rois','roiNames');

% bar( shTmaps);title('t split half corr (across vox)','interpreter','Latex'); hold all;
% h = gca; set(h,'FontSize',12);set(h,'XTickLabel',connamesNoUnderscore,'TickLabelInterpreter','Latex');%

%% 