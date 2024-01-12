function [PP_Bn, PP_An, PP_AnNS] = PairWiseCorrelation(tonoBn, tonoAn, tonoAnNS, roi_large, roi, roiID, folderOut, subj)

PPCor_Bn = [];
PPCor_An = [];
PPCor_AnNS = [];
PP_Bn = [];
PP_An = [];
PP_AnNS = [];

uroi = unique(roiID);

for i = 1:size(uroi,2)

    index = roi(find(roiID==i));

    PPCor_Bn = corr(tonoBn(index,:));
    PPCor_An = corr(tonoAn(index,:));
    PPCor_AnNS = corr(tonoAnNS(index,:));
    PPCor_Bn = triu(PPCor_Bn,1);
    PPCor_An = triu(PPCor_An,1);
    PPCor_AnNS = triu(PPCor_AnNS,1);
    PP_Bn(:,i) = PPCor_Bn(find(PPCor_Bn~=0));
    PP_An(:,i) = PPCor_An(find(PPCor_An~=0));
    PP_AnNS(:,i) = PPCor_AnNS(find(PPCor_AnNS~=0));

    roiCount(i,:,:) = repmat(i,[1,size(PP_Bn,1),3]);
end 

methCount(:,:,1) = repmat(1, [5, size(PP_Bn,1)]);
methCount(:,:,2) = repmat(2, [5, size(PP_Bn,1)]);
methCount(:,:,3) = repmat(3, [5, size(PP_Bn,1)]);

PP_Bn = PP_Bn';
PP_An = PP_An';
PP_AnNS = PP_AnNS';

data(:,:,1) = PP_Bn;
data(:,:,2) = PP_An;
data(:,:,3) = PP_AnNS;

figure
boxchart(2.*roiCount(:), data(:), 'GroupByColor', methCount(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none')
title('Pairwise correlation per ROI')
xlabel('ROIs')
ylabel('Correlation value')
set(gca, 'XTick', [])

s = gcf;
saveas(s,[folderOut, [subj,'_PairWiseCorr.fig']]);
save([folderOut,[subj, '_PairWiseCorr.mat']],'data');
saveas(s, [folderOut, [subj,'_PairWiseCorr.svg']]);

