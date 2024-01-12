function [tonoBN, tonoAN] = get_tono_wholeBstd(colnames, tonoconH,tonoconL,tbyCond,tbyCondSH1,tbyCondSH2,bBn,bAn,splithalf_exp,roiID,ROIname,p0,trial_x_run)
disp(['Computing tonotopies']);
% whole brain standardization
% tonoconH = {'Comp_H','Oddball_H' ,'Unexp_H'};
% tonoconL = {'Comp_L','Oddball_L' ,'Unexp_L'};
mtalltrials    = squeeze(mean(tbyCond(1:2,:,:),1));% for selection we use only Comp_H and Comp_L. tbyCond is 6 cond x 2 (before and after Nordic) x Nvox
mtalltrialsSH1 = squeeze(mean(tbyCondSH1(1:2,:,:),1));% for selection we use only Comp_H and Comp_L. tbyCond is 6 cond x 2 (before and after Nordic) x Nvox
mtalltrialsSH2 = squeeze(mean(tbyCondSH2(1:2,:,:),1));% for selection we use only Comp_H and Comp_L. tbyCond is 6 cond x 2 (before and after Nordic) x Nvox
trial_x_run_dummy = dummyvar(trial_x_run);
% p0         = 90;
nvoxROI    = size(bBn,2);
% unique integers defining ROIs
uroi      = unique(roiID); nroi = length(uroi);
nrep      = 10;
nruns     = size(trial_x_run_dummy,2);
for itrep = 1:nrep
    runsperm = randperm(nruns);
    splithalf_exp  = sum(trial_x_run_dummy(:,runsperm(1:nruns/2)),2)';
    for itroi = 1:nroi%  loop across rois
        roisLoc = roiID == uroi(itroi) ; % a binary vector indicating the precense of ROI
        nvoxROI(itroi)   = sum(roisLoc);
        for itcon = 1:length(tonoconL)
            selBetasH = contains(colnames , tonoconH{itcon});
            selBetasL = contains(colnames , tonoconL{itcon});
            [tonoBN{itroi}(itcon,:) ,   tonoAN{itroi}(itcon,:) ,   tonoBNsel{itroi}(itcon,:) ,   tonoANsel{itroi}(itcon,:)     ] = sub_tono_before_after(     selBetasL,  selBetasH,  bBn,bAn, 1,roisLoc,mtalltrials,p0 );
            [tonoBNsh1{itroi}(itcon,itrep,:) ,tonoANsh1{itroi}(itcon,itrep,:) ,tonoBNsh1SEL{itroi}(itcon,itrep,:) ,tonoANsh1SEL{itroi}(itcon,itrep,:)  ] = sub_tono_before_after( splithalf_exp & selBetasL,   splithalf_exp & selBetasH, bBn,bAn, 1,roisLoc,mtalltrials,p0 );
            [tonoBNsh2{itroi}(itcon,itrep,:) ,tonoANsh2{itroi}(itcon,itrep,:) ,tonoBNsh2SEL{itroi}(itcon,itrep,:) ,tonoANsh2SEL{itroi}(itcon,itrep,:)  ] = sub_tono_before_after(~splithalf_exp & selBetasL, ~splithalf_exp & selBetasH,  bBn,bAn, 1,roisLoc,mtalltrials,p0 );
        end
    end
end
%%
methID = repmat([1 1 1 2 2 2]',1,nrep);condID = repmat([1 2 3 1 2 3 ]',1,nrep);
figure('Color',[1,1,1],'Position',[281 19 1160 820]),
for itroi = 1:nroi
    
    subplot(nroi, 2, 2*(itroi-1)+1);
    tonoOverlap    = [sum(tonoBNsh1{itroi} ==  tonoBNsh2{itroi},3)/nvoxROI(itroi)  ;  sum(tonoANsh1{itroi} ==  tonoANsh2{itroi},3)/nvoxROI(itroi) ];
    tonoOverlap    = squeeze(tonoOverlap);
    boxchart(condID(:),tonoOverlap(:),'GroupByColor',methID(:))
    ylim([0.5 0.6]);hold all; %xlim([0.8 3.2])
    low  = mean(sum( tonoBNsh1{itroi}== 1 ,3)/nvoxROI(itroi) ,2);
    high = mean(sum( tonoBNsh1{itroi}==-1,3)/nvoxROI(itroi)  ,2);
    plot( xlim(), [max(high,low) max(high,low)],'k--' );
    title([ROIname{itroi}, ' no vox selection'],'interpreter','Latex')
    h = gca; set(h,'FontSize',12);
    set(h,'XTick',condID(1:3,1))
    set(h,'XTickLabel',{'Comp','Oddball','Unexp'},'TickLabelInterpreter','Latex')
    
    subplot(nroi,2,2*(itroi-1)+2)
    tonoOverlapSel = [ sum(tonoBNsh1SEL{itroi} ==  tonoBNsh2SEL{itroi}, 3) /size(tonoBNsh1SEL{itroi},3)  ;  sum(tonoANsh1SEL{itroi} ==  tonoANsh2SEL{itroi}, 3) /size(tonoANsh1SEL{itroi},3)   ];
    tonoOverlapSel = squeeze(tonoOverlapSel);
%    boxplot( tonoOverlapSel');
%    bar(tonoOverlapSel');
%    bar([mean(tonoOverlap(1:3,:),2)  mean(tonoOverlap(4:6,2),2)])
     boxchart(condID(:),tonoOverlapSel(:),'GroupByColor',methID(:))
    ylim([0.5 0.6]);hold all
    title([ROIname{itroi},' selection with F p0 = ', num2str(p0) ],'interpreter','Latex')
    hold all;
    low  = mean( sum( tonoBNsh1SEL{itroi} == 1 ,3)/size(tonoBNsh1SEL{itroi},3),2);
    high = mean( sum( tonoBNsh1SEL{itroi} ==-1 ,3)/size(tonoBNsh1SEL{itroi},3),2);
    plot( xlim(), [max(high,low) max(high,low)],'k--' );
%     ylim([0.45 0.6]);
    h = gca; set(h,'FontSize',12);
     set(h,'XTick',condID(1:3,1))
    set(h,'XTickLabel',{'Comp','Oddball','Unexp'},'TickLabelInterpreter','Latex')
end
