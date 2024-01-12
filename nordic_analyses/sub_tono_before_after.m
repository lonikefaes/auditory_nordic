function [tonoBN,tonoAN,tonoBNsel,tonoANsel] = sub_tono_before_after(selBetasL, selBetasH, bBn,bAn, globalZscore,roisLoc,tmap,p0 )
% thre
tTH       = prctile(tmap',p0);
% average x cond
        mHbn      = (mean(bBn(selBetasH,:))); % mean across trials of the corresponding condition
        mLbn      = (mean(bBn(selBetasL,:))); % standardization force each trial to have zero mean and variance 1 across the whole ROI (brain)
        mHan      = (mean(bAn(selBetasH,:))); % mean across trials of the corresponding condition
        mLan      = (mean(bAn(selBetasL,:)));
%std
if(globalZscore)
        cortexmeanHbn = mean(mHbn);  cortexvarHbn  = var(mHbn);  %High before mean and variance without looking at the ROI (whole cortex) 
        cortexmeanLbn = mean(mLbn);  cortexvarLbn  = var(mLbn);  %Low before mean and variance without looking at the ROI (whole cortex) 
        cortexmeanHan = mean(mHan);  cortexvarHan  = var(mHan);  %High before mean and variance without looking at the ROI (whole cortex) 
        cortexmeanLan = mean(mLan);  cortexvarLan  = var(mLan);  %Low before mean and variance without looking at the ROI (whole cortex) 


        cortexmeanHbnSel = mean(mHbn(tmap(1,:) >  tTH(1) ));  cortexvarHbnSel  = var(mHbn(tmap(1,:) >  tTH(1) ));  %High before mean and variance without looking at the ROI (whole cortex) 
        cortexmeanLbnSel = mean(mLbn(tmap(1,:) >  tTH(1) ));  cortexvarLbnSel  = var(mLbn(tmap(1,:) >  tTH(1) ));  %Low before mean and variance without looking at the ROI (whole cortex) 
        cortexmeanHanSel = mean(mHan(tmap(2,:) >  tTH(2) ));  cortexvarHanSel  = var(mHan(tmap(2,:) >  tTH(2) ));  %High before mean and variance without looking at the ROI (whole cortex) 
        cortexmeanLanSel = mean(mLan(tmap(2,:) >  tTH(2) ));  cortexvarLanSel  = var(mLan(tmap(2,:) >  tTH(2) ));  %Low before mean and variance without looking at the ROI (whole cortex) 
end
% tono no selection
        tonoBN = sign( (mLbn(roisLoc) - cortexmeanLbn)./sqrt(cortexvarLbn)  - (mHbn(roisLoc) - cortexmeanHbn)./sqrt(cortexvarHbn) ); % Low - High
        tonoAN = sign( (mLan(roisLoc) - cortexmeanLan)./sqrt(cortexvarLan)  - (mHan(roisLoc) - cortexmeanHan)./sqrt(cortexvarHan) ); % Low - High
 % tono with selection        
        tonoBNsel = sign( (mLbn(roisLoc & tmap(1,:) >  tTH(1)) - cortexmeanLbnSel)./sqrt(cortexvarLbnSel)  - (mHbn(roisLoc & tmap(1,:) >  tTH(1))  - cortexmeanHbnSel)./sqrt(cortexvarHbnSel) ); % Low - High
        tonoANsel = sign( (mLan(roisLoc & tmap(2,:) >  tTH(2)) - cortexmeanLanSel)./sqrt(cortexvarLanSel)  - (mHan(roisLoc & tmap(2,:) >  tTH(2))  - cortexmeanHanSel)./sqrt(cortexvarHanSel) ); % Low - High        
