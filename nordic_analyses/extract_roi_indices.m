function [roi , roiID] = extract_roi_indices(roinames,datafolder,subj)
% extract roi indices and identificators integers from files
roi       = [];roiID     = [];
for itroi = 1:length(roinames)
    targetroi   = roinames{itroi};%'HG_GM.msk';% set the ROI name here
     roiT = xff([datafolder,subj,filesep,'ROI',filesep,targetroi]);
%       roiT = xff([datafolder,filesep,filesep,targetroi]);

    roiT = roiT.Mask; roiT = find(roiT);
    roi  = [roi roiT'];
    roiID = [roiID  itroi + 0.*roiT'];
    disp(['loading ', num2str(roinames{itroi}), ' with nvox ',num2str(length(roiT))]);
end

