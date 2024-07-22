function [roiOut , roiIDOut] = IntersectROIs(roi,roiID,roi_large)

%%% this function takes the indices of a large ROI and small ROIs and
%%% intersects them

nrROIs = max(roiID);
roiOut = [];
roiIDOut = [];
for i=1:nrROIs
    this_roi = roi(find(roiID==i));
    [~,index,~] = intersect (roi_large,this_roi);
    roiOut = [roiOut,index'];
    roiIDOut = [roiIDOut,i*ones(1,length(index))];
end