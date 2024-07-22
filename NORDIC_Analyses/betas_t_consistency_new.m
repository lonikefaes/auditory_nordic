    function [tbyCond,tbyCondSH1,tbyCondSH2,BetasCond ] = betas_t_consistency_new(colnames,connames,connamesNoUnderscore,bBn,bAn,bAnNS,splithalf_exp , roiID, trial_x_run, nrep,roiNames,p0,folderOut,subj)
    % Do the reliability analysis of betas and t 
    % bBn betas before NORDIC
    % bAn betas after NORDIC
    % connames             = {'Comp_H' , 'Comp_L', 'Oddball_L','Oddball_H', 'Unexp_H', 'Unexp_L' };
    % connamesNoUnderscore = {'Comp H' , 'Comp L', 'Oddball L','Oddball H', 'Unexp H', 'Unexp L' };
    % splithalf_exp     = zeros(size(colnames));% this vector defines which trials are in the first half of the experiment
    % splithalf_exp( 1 : length(splithalf_exp)/2)=1;
    
        for itcon = 1:length(connames)  
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
             BetasCondSH1(itcon,:)   =  mean( bBn(con & splithalf_exp,:) ); % Taking the mean betas from the trials that are divided in two halfs
             BetasCondSH2(itcon,:)   =  mean( bBn(con & ~splithalf_exp,:) );   
             
        end 
    
        %These things are not used atm
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
            end
         relBetas(itroi,itrep,:) = [ corr(mbBsh1(locroi)' ,mbBsh2(locroi)' )  corr(mbAsh1(locroi)' ,mbAsh2(locroi)' )  corr(mbAnssh1(locroi)' ,mbAnssh2(locroi)' )  ];
         relTs(itroi,itrep,:)    = [ corr(STATSbnSH1.tstat(locroi)' ,STATSbnSH2.tstat(locroi)' )  corr(STATSanSH1.tstat(locroi)' ,STATSanSH2.tstat(locroi)' )  corr(STATSanNSSH1.tstat(locroi)' ,STATSanNSSH2.tstat(locroi)' )];
         roiCount(itroi,itrep,:) = [itroi,itroi,itroi] ;
         methCount(itroi,itrep,:)= [1 2 3];
         BetasCondrepSH(itroi,itrep,:) =  [median( [  mbBsh2(locroi)' ])      median( mbAsh2(locroi)')       median( mbAnssh2(locroi)')  ]; %We take the median % beta change of one of the splits and repeat 50 times.
         tCondrepSH(itroi,itrep,:)     =  [median(  STATSbnSH2.tstat(locroi)')           median(STATSanSH2.tstat(locroi)')   median(STATSanNSSH2.tstat(locroi)')   ];
        end
    end

    
    %% Original figure of beta change, t and reliability
    disp(['Creating the original figure of beta and t change and realibility' ]);



    rois = 2.*uroi;
    reliability_and_tmaps(rois,roiCount,roiNames,methCount,BetasCondrepSH, tCondrepSH , relBetas, relTs); % this is just the code coming next pack in a function
%     figure('Color',[1,1,1],'Position',[81 59 2060 820]),
%     subplot(221)
%     boxchart(2.*roiCount(:),BetasCondrepSH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
%     
%     % mb  = bar(mBetasCond); 
%     % legend('before','after'); title('$\beta$','interpreter','Latex');
%     % ylim([prctile(mBetasCompleteLH(:),10),  prctile(mBetasCompleteLH(:),90)])
%     ylim([-1.5 3]);
%     h   = gca; set(h,'FontSize',12); 
%     set(h,'XTick',unique(rois(:)))
%     set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
%     ylabel('$\beta$  \% change','interpreter','Latex')
%     
%     subplot(222);
%     mtCompleteLH = squeeze(mean(tbyCond(1:2,:,:),1));
%     boxchart(2.*roiCount(:),tCondrepSH(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
%     hold all;
%     h = gca; set(h,'FontSize',12);
%     % ylim([prctile(tCondrepSH(:),1), prctile(tCondrepSH(:),99)])
%     ylim([-1 3]);
%     set(h,'XTick',unique(rois(:)))
%     set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
%     ylabel('t-statistic','interpreter','Latex')
%     legend('Original', 'NORDIC', 'NORDIC No Noise', 'Location', 'northeast')
% 
%     
%     % errorbar(mb(1,1).XEndPoints,mtbyCond(:,1),sdCond(:,1),'LineWidth',2)
%     subplot(223);
%     boxchart(2.*roiCount(:),relBetas(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
%     
%     title('$\beta$ reliability','interpreter','Latex'); hold all;
%     h = gca; set(h,'FontSize',12);
%     ylim([-0.2 1])
%     set(h,'XTick',unique(rois(:)))
%     set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
%     ylabel('split half corr','interpreter','Latex')
%     % bar( shBetas);title('$\beta$ split half corr (across vox)','interpreter','Latex'); hold all;
%     % h = gca; set(h,'FontSize',12);set(h,'XTickLabel',connamesNoUnderscore,'TickLabelInterpreter','Latex')
%     
%     subplot(224);
%     boxchart(2.*roiCount(:),relTs(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
%     
%     title('t-value reliability','interpreter','Latex'); hold all;
%     h = gca; set(h,'FontSize',12);
%     ylim([-0.2 1])
%     set(h,'XTick',unique(rois(:)))
%     set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
%     ylabel('split half corr','interpreter','Latex')
%     
    h = gcf;

    saveas(h,[folderOut, [subj,'_BetasANDt_changeANDReliability.fig']]);
    saveas(h, [folderOut, [subj,'_BetasANDt_changeANDReliability.svg']]);
    save([folderOut,[subj, '_BetasANDt_changeANDReliability.mat']],'uroi','roiCount','BetasCondrepSH','methCount','relBetas','relTs','rois','roiNames', 'tCondrepSH');
    
    % bar( shTmaps);title('t split half corr (across vox)','interpreter','Latex'); hold all;
    % h = gca; set(h,'FontSize',12);set(h,'XTickLabel',connamesNoUnderscore,'TickLabelInterpreter','Latex');%
    
   
    %% Calculating betas and t-value change per condition
    disp(['Calculating split halfs' ]);


    uroi = unique(roiID);
    trial_x_run_dummy = dummyvar(trial_x_run);
    nruns             = size(trial_x_run_dummy,2);

    relBetas_perCond_Bn = zeros(length(uroi), nrep, nruns, 3);
    relBetas_perCond_An = zeros(length(uroi), nrep, nruns, 3);
    relBetas_perCond_AnNs = zeros(length(uroi), nrep, nruns, 3);
    BetasCondrepSH_perCond_Bn = zeros(length(uroi), nrep, nruns, 3);
    BetasCondrepSH_perCond_An = zeros(length(uroi), nrep, nruns, 3);
    BetasCondrepSH_perCond_AnNs = zeros(length(uroi), nrep, nruns, 3);
    tCondrepSH_perCond_Bn = zeros(length(uroi), nrep, nruns, 3);
    tCondrepSH_perCond_An = zeros(length(uroi), nrep, nruns, 3);
    tCondrepSH_perCond_AnNs = zeros(length(uroi), nrep, nruns, 3);

    for itrep = 1:nrep
        runsperm       = randperm(nruns);
        splithalf_exp  = sum(trial_x_run_dummy(:,runsperm(1:nruns/2)),2)';


        for itcon = 1:length(connames)
            cons = contains(colnames, connames{itcon});% All conditions
            %     %no split variability
            mbBsh1_perCond = mean(bBn(cons==1 &  splithalf_exp,:));
            mbBsh2_perCond = mean(bBn(cons==1 & ~splithalf_exp,:));
            mbAsh1_perCond = mean(bAn(cons==1 &  splithalf_exp,:));
            mbAsh2_perCond = mean(bAn(cons==1 & ~splithalf_exp,:));
            mbAnssh1_perCond = mean(bAnNS(cons==1 &  splithalf_exp,:));
            mbAnssh2_perCond = mean(bAnNS(cons==1 & ~splithalf_exp,:));

            [H,P,CI,STATSbnSH1_perCond]    = ttest(bBn(cons==1 & splithalf_exp,:)  );     [H,P,CI,STATSbnSH2_perCond] = ttest(bBn(cons==1 & ~splithalf_exp,:)  );
            [H,P,CI,STATSanSH1_perCond]    = ttest(bAn(cons==1 & splithalf_exp,:)  );     [H,P,CI,STATSanSH2_perCond] = ttest(bAn(cons==1 & ~splithalf_exp,:)  );
            [H,P,CI,STATSanNSSH1_perCond]  = ttest(bAnNS(cons==1 & splithalf_exp,:)  );     [H,P,CI,STATSanNSSH2_perCond] = ttest(bAnNS(cons==1 & ~splithalf_exp,:)  );


            for itroi = 1:length(uroi)
                locroi = find(roiID  == uroi(itroi));
                if(p0 > 0)  % select those voxels above a percentile thres at each ROI
                    selvox = STATSbnSH1_perCond.tstat( locroi ) > prctile( STATSbnSH1_perCond.tstat( locroi ) , p0);
                    locroi = locroi(selvox);
                end

                relBetas_perCond_Bn(itroi,itrep,itcon,:) = [corr(mbBsh1_perCond(locroi)' ,mbBsh2_perCond(locroi)') corr(mbAsh1_perCond(locroi)' ,mbAsh2_perCond(locroi)') corr(mbAnssh1_perCond(locroi)' ,mbAnssh2_perCond(locroi)')];
                roiCount_perCond(itroi,itrep,:,:) = [itroi, itroi, itroi, itroi, itroi, itroi] ;
                condCount(itroi,itrep,:) = [1 2 3 4 5 6];
                BetasCondrepSH_perCond_Bn(itroi,itrep,itcon,:) = [median(mbBsh2_perCond(locroi)') median(mbAsh2_perCond(locroi)') median(mbAnssh2_perCond(locroi)')];
                %BetasCondrepSH_perCond_An(itroi,itrep,itcon,:) = median(mbAsh2_perCond(locroi)');
                %BetasCondrepSH_perCond_AnNS(itroi,itrep,itcon,:) = median(mbAnssh2_perCond(locroi)');

                tCondrepSH_perCond_Bn(itroi,itrep,itcon,:) =  [median(STATSbnSH2_perCond.tstat(locroi)') median(STATSanSH2_perCond.tstat(locroi)') median(STATSanNSSH2_perCond.tstat(locroi)')];
                %tCondrepSH_perCond_An(itroi,itrep,itcon) = median(STATSanSH2_perCond.tstat(locroi)');
                %tCondrepSH_perCond_AnNS(itroi,itrep,itcon) = median(STATSanNSSH2_perCond.tstat(locroi)');

            end
        end
    end


    %% Plot the betas per condition and ROI

    disp(['Plotting beta change per condition and per ROI' ]);


    BetasCompH = squeeze(BetasCondrepSH_perCond_Bn(:,:,1,:));
    BetasCompL = squeeze(BetasCondrepSH_perCond_Bn(:,:,2,:));
    BetasOddL = squeeze(BetasCondrepSH_perCond_Bn(:,:,3,:));
    BetasOddH = squeeze(BetasCondrepSH_perCond_Bn(:,:,4,:));
    BetasUnexpH = squeeze(BetasCondrepSH_perCond_Bn(:,:,5,:));
    BetasUnexpL = squeeze(BetasCondrepSH_perCond_Bn(:,:,6,:));


    figure('Color',[1,1,1],'Position',[81 59 2060 820]),
    sgtitle('Beta change per condition and ROI', 'interpreter','Latex')
    subplot(321)
    boxchart(2.*roiCount(:),BetasCompH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 3]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    title('PredH','interpreter','Latex')
    ylabel('$\beta$  \% change','interpreter','Latex')

    subplot(322)
    boxchart(2.*roiCount(:),BetasCompL(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 3]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('$\beta$  \% change','interpreter','Latex')
    title('PredL','interpreter','Latex')
    legend('Original', 'NORDIC', 'NORDIC No Noise', 'Location', 'northeast')


    subplot(323)
    boxchart(2.*roiCount(:),BetasOddH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 3]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('$\beta$  \% change','interpreter','Latex')
    title('PredH','interpreter','Latex')

    subplot(324)
    boxchart(2.*roiCount(:),BetasOddL(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 3]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('$\beta$  \% change','interpreter','Latex')
    title('PredL','interpreter','Latex')

    subplot(325)
    boxchart(2.*roiCount(:),BetasUnexpH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 3]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('$\beta$  \% change','interpreter','Latex')
    title('UnpredH','interpreter','Latex')

    subplot(326)
    boxchart(2.*roiCount(:),BetasUnexpL(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 3]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('$\beta$  \% change','interpreter','Latex')
    title('UnpredL','interpreter','Latex')

    h = gcf;
    saveas(h,[folderOut, [subj,'_BetasperCond_ROI.fig']]);
    saveas(h, [folderOut, [subj,'_BetasperCond_ROI.svg']]);
    save([folderOut,[subj, '_BetasperCond_ROI.mat']],'uroi','roiCount','BetasCompH','BetasCompL','BetasOddH','BetasOddL','BetasUnexpH','BetasUnexpL', 'methCount', 'roiNames');

%% Plot t-values for each condition
    
    disp(['Plotting t change per condition and per ROI' ]);


    tCompH = squeeze(tCondrepSH_perCond_Bn(:,:,1,:));
    tCompL = squeeze(tCondrepSH_perCond_Bn(:,:,2,:));
    tOddL = squeeze(tCondrepSH_perCond_Bn(:,:,3,:));
    tOddH = squeeze(tCondrepSH_perCond_Bn(:,:,4,:));
    tUnexpH = squeeze(tCondrepSH_perCond_Bn(:,:,5,:));
    tUnexpL = squeeze(tCondrepSH_perCond_Bn(:,:,6,:));


    figure('Color',[1,1,1],'Position',[81 59 2060 820]), hold on,
    sgtitle('t value change per condition and ROI', 'interpreter','Latex')

    subplot(321)
    boxchart(2.*roiCount(:),tCompH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 2]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    title('PredH','interpreter','Latex')
    ylabel('t-value change','interpreter','Latex')

    subplot(322)
    boxchart(2.*roiCount(:),tCompL(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 2]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('t-value change','interpreter','Latex')
    title('PredL','interpreter','Latex')
    legend('Original', 'NORDIC', 'NORDIC No Noise', 'Location', 'northeast')


    subplot(323)
    boxchart(2.*roiCount(:),tOddH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 2]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('t-value change','interpreter','Latex')
    title('MispredH','interpreter','Latex')

    subplot(324)
    boxchart(2.*roiCount(:),tOddL(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 2]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('t-value change','interpreter','Latex')
    title('MispredL','interpreter','Latex')

    subplot(325)
    boxchart(2.*roiCount(:),tUnexpH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 2]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('t-value change','interpreter','Latex')
    title('UnpredH','interpreter','Latex')

    subplot(326)
    boxchart(2.*roiCount(:),tUnexpL(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    ylim([0 2]);

    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('t-value change','interpreter','Latex')
    title('UnpredL','interpreter','Latex')


    h = gcf;

    saveas(h,[folderOut, [subj,'_TperCond_ROI.fig']]);
    saveas(h, [folderOut, [subj,'_TperCond_ROI.svg']]);
    save([folderOut,[subj, '_TperCond_ROI.mat']],'uroi','roiCount','tCompH','tCompL','tOddH','tOddL','tUnexpH','tUnexpL', 'methCount', 'roiNames');

