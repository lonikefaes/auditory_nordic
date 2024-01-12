function [] = reliability_and_tmaps(rois,roiCount,roiNames,methCount,BetasCondrepSH, tCondrepSH , relBetas, relTs)

%     rois = 2.*uroi;
    figure('Color',[1,1,1], 'Position',[281 19 1000 700])
    subplot(221)
    mb = boxchart(2.*roiCount(:),BetasCondrepSH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    % mb  = bar(mBetasCond); 
    % legend('before','after'); title('$\beta$','interpreter','Latex');
    % ylim([prctile(mBetasCompleteLH(:),10),  prctile(mBetasCompleteLH(:),90)])
    ylim([-1 4]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('$\beta$  \% change','interpreter','Latex')


    subplot(222);
%     mtCompleteLH = squeeze(mean(tbyCond(1:2,:,:),1));
    boxchart(2.*roiCount(:),tCondrepSH(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
    hold all;
    h = gca; set(h,'FontSize',12);
    ylim([prctile(tCondrepSH(:),1), prctile(tCondrepSH(:),99)])
    ylim([-1 8]);
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('t-statistic','interpreter','Latex')
    legend('Original', 'NORdef', 'NORnn', 'Location', 'northeast', 'FontSize', 10, 'FontName', 'Arial')


    
    % errorbar(mb(1,1).XEndPoints,mtbyCond(:,1),sdCond(:,1),'LineWidth',2)
    subplot(223);
    boxchart(2.*roiCount(:),relBetas(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
    
    title('$\beta$ reliability','interpreter','Latex'); hold all;
    h = gca; set(h,'FontSize',12);
    ylim([-0.2 1])
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('split half corr','interpreter','Latex')
  

    
    subplot(224);
    boxchart(2.*roiCount(:),relTs(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
    
    title('t-value reliability','interpreter','Latex'); hold all;
    h = gca; set(h,'FontSize',12);
    ylim([-0.2 1])
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('split half corr','interpreter','Latex')
