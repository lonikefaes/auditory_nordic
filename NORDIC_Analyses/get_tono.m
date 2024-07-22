function [tonoBN, tonoAN,rsaBNAN  ,rsaBNANnn  ,rsaANANnn,prefBN,prefAN,prefANnn] = get_tono(colnames, tonoconH,tonoconL,tbyCond,bBn,bAn,bAnNS,splithalf_exp,roiID,ROIname,trial_x_run,p0, folderOut)
disp(['Computing tonotopies']);
% within ROI standardization
% tonoconH = {'Comp_H','Oddball_H' ,'Unexp_H'};
% tonoconL = {'Comp_L','Oddball_L' ,'Unexp_L'};
mtalltrials = squeeze(mean(tbyCond(1:2,:,:),1));% for selection we use only Comp_H and Comp_L
% p0        = 90;
nvoxROI   = size(bBn,2);
% unique integers defining ROIs
uroi      = unique(roiID); nroi = length(uroi);
%%
trial_x_run_dummy = dummyvar(trial_x_run);
nruns     = size(trial_x_run_dummy,2);
nrep      = 30;
for itrep = 1:nrep
    runsperm       = randperm(nruns);
    splithalf_exp  = sum( trial_x_run_dummy(:, runsperm(1:nruns/2) ),2)';

    for itroi = 1:nroi%  loop across rois
        roisLoc = roiID == uroi(itroi) ; % a binary vector indicating the presence of ROI
        nvoxROI(itroi)   = sum(roisLoc);
        tTH       = prctile(mtalltrials(:,roisLoc)',p0);% gets the t value that makes the p0 perctentile in this ROIs
        
        for itcon = 1:length(tonoconL)
            selBetasH = contains(colnames , tonoconH{itcon});
            selBetasL = contains(colnames , tonoconL{itcon});
            mHbn      = (mean(bBn(selBetasH,:))); % mean across trials of the corresponding condition
            mLbn      = (mean(bBn(selBetasL,:))); % standardization force each trial to have zero mean and variance 1 across the whole ROI (brain)
            mHan      = (mean(bAn(selBetasH,:))); % mean across trials of the corresponding condition
            mLan      = (mean(bAn(selBetasL,:)));
            mLanNN    = (mean(bAnNS(selBetasL,:)));
            mHanNN    = (mean(bAnNS(selBetasH,:)));
            % zscoring within ROI with no selection
            tonoBN{itcon,itroi}   = (zscore(mLbn(mtalltrials(1,:)>tTH(1) & roisLoc))   - zscore(mHbn(mtalltrials(1,:)>tTH(1) & roisLoc)));%./(zscore(mLbn(mtalltrials(1,:)>tTH(1) & roisLoc))   + zscore(mHbn(mtalltrials(1,:)>tTH(1) & roisLoc))); % Low - High
            tonoAN{itcon,itroi}   = (zscore(mLan(mtalltrials(1,:)>tTH(1) & roisLoc))   - zscore(mHan(mtalltrials(1,:)>tTH(1) & roisLoc)));%./(zscore(mLan(mtalltrials(2,:)>tTH(2) & roisLoc))   + zscore(mHan(mtalltrials(2,:)>tTH(2) & roisLoc))); % Low - High
            tonoANnn{itcon,itroi} = (zscore(mLanNN(mtalltrials(1,:)>tTH(1) & roisLoc)) - zscore(mHanNN(mtalltrials(1,:)>tTH(1) & roisLoc)));%./(zscore(mLanNN(mtalltrials(3,:)>tTH(3) & roisLoc)) + zscore(mHanNN(mtalltrials(3,:)>tTH(3) & roisLoc))); % Low - High

%             tonoBN{itcon,itroi}   = tonoangle( mLbn(mtalltrials(1,:)>tTH(1) & roisLoc) ,   mHbn(mtalltrials(1,:)>tTH(1) & roisLoc) ); % Low - High

            % split half tonotopy maps: first split
            mHbn1      = (mean(bBn(splithalf_exp & selBetasH,:))); % mean across trials of the corresponding condition
            mLbn1      = (mean(bBn(splithalf_exp & selBetasL,:))); % standardization force each trial to have zero mean and variance 1 across the whole ROI (brain)
            mHan1      = (mean(bAn(splithalf_exp & selBetasH,:))); % mean across trials of the corresponding condition
            mLan1      = (mean(bAn(splithalf_exp & selBetasL,:)));
            mHanNS1    = (mean(bAnNS(splithalf_exp & selBetasH,:))); % mean across trials of the corresponding condition
            mLanNS1    = (mean(bAnNS(splithalf_exp & selBetasL,:)));
       
            % split half tonotopy maps: second split
            mHbn2      = (mean(bBn(~splithalf_exp & selBetasH,:))); % mean across trials of the corresponding condition
            mLbn2      = (mean(bBn(~splithalf_exp & selBetasL,:))); % standardization force each trial to have zero mean and variance 1 across the whole ROI (brain)
            mHan2      = (mean(bAn(~splithalf_exp & selBetasH,:))); % mean across trials of the corresponding condition
            mLan2      = (mean(bAn(~splithalf_exp & selBetasL,:)));
            mHanNS2    = (mean(bAnNS(~splithalf_exp & selBetasH,:))); % mean across trials of the corresponding condition
            mLanNS2    = (mean(bAnNS(~splithalf_exp & selBetasL,:)));
 
            tonoBNsh1{itroi}(itcon,itrep,:)   = sign(zscore(mLbn1(roisLoc) )   -  zscore(mHbn1(roisLoc))); % Low - High
            tonoANsh1{itroi}(itcon,itrep,:)   = sign(zscore(mLan1(roisLoc))    -  zscore(mHan1(roisLoc))); % Low - High
            tonoANnssh1{itroi}(itcon,itrep,:) = sign(zscore(mLanNS1(roisLoc))  -  zscore(mHanNS1(roisLoc))); % Low - High
     
            tonoBNsh2{itroi}(itcon,itrep,:)   = sign(zscore(mLbn2(roisLoc))    -  zscore(mHbn2(roisLoc))); % Low - High
            tonoANsh2{itroi}(itcon,itrep,:)   = sign(zscore(mLan2(roisLoc))    -  zscore(mHan2(roisLoc))); % Low - High
            tonoANnssh2{itroi}(itcon,itrep,:) = sign(zscore(mLanNS2(roisLoc))  -  zscore(mHanNS2(roisLoc))); % Low - High            
            
            tonoBNsh1SEL{itroi}(itcon,itrep,:)   = sign( zscore(mLbn1( mtalltrials(1,:)>tTH(1) & roisLoc )) - zscore(mHbn1( mtalltrials(1,:)>tTH(1) & roisLoc ))  ); % Low - High
            tonoANsh1SEL{itroi}(itcon,itrep,:)   = sign( zscore(mLan1( mtalltrials(2,:)>tTH(2) & roisLoc )) - zscore(mHan1( mtalltrials(2,:)>tTH(2) & roisLoc ))  ); % Low - High
            tonoANnssh1SEL{itroi}(itcon,itrep,:) = sign( zscore(mLanNS1( mtalltrials(3,:)>tTH(3) & roisLoc )) - zscore(mHanNS1( mtalltrials(3,:)>tTH(3) & roisLoc ))  ); % Low - High

            tonoBNsh2SEL{itroi}(itcon,itrep,:)   = sign( zscore(mLbn2( mtalltrials(1,:)>tTH(1) & roisLoc )) - zscore(mHbn2( mtalltrials(1,:)>tTH(1) & roisLoc ))  ); % Low - High
            tonoANsh2SEL{itroi}(itcon,itrep,:)   = sign( zscore(mLan2( mtalltrials(2,:)>tTH(2) & roisLoc )) - zscore(mHan2( mtalltrials(2,:)>tTH(2) & roisLoc ))  ); % Low - High
            tonoANnssh2SEL{itroi}(itcon,itrep,:) = sign( zscore(mLanNS2( mtalltrials(3,:)>tTH(3) & roisLoc )) - zscore(mHanNS2( mtalltrials(3,:)>tTH(3) & roisLoc ))  ); % Low - High
            
            
            % case circular tonotopy
%             tonoBNsh1{itroi}(itcon,itrep,:)     = tonoangle(mLbn1(roisLoc)   , mHbn1(roisLoc)); % Low - High
%             tonoANsh1{itroi}(itcon,itrep,:)     = tonoangle(mLan1(roisLoc)   , mHan1(roisLoc)); % Low - High
%             tonoANnssh1{itroi}(itcon,itrep,:)   = tonoangle(mLanNS1(roisLoc)   , mHanNS1(roisLoc)); % Low - High
%      
%             tonoBNsh2{itroi}(itcon,itrep,:)      =  tonoangle(mLbn2(roisLoc)     , mHbn2(roisLoc)); % Low - High
%             tonoANsh2{itroi}(itcon,itrep,:)      =  tonoangle(mLan2(roisLoc)     , mHan2(roisLoc));% Low - High
%             tonoANnssh2{itroi}(itcon,itrep,:)    =  tonoangle(mLanNS2(roisLoc)   , mHanNS2(roisLoc)); % Low - High            
%             
%             tonoBNsh1SEL{itroi}(itcon,itrep,:)   =  tonoangle(mLbn1( mtalltrials(1,:)>tTH(1) & roisLoc )   , mHbn1( mtalltrials(1,:)>tTH(1) & roisLoc )  ); % Low - High
%             tonoANsh1SEL{itroi}(itcon,itrep,:)   =  tonoangle(mLan1( mtalltrials(2,:)>tTH(2) & roisLoc )   , mHan1( mtalltrials(2,:)>tTH(2) & roisLoc )  ); % Low - High
%             tonoANnssh1SEL{itroi}(itcon,itrep,:) =  tonoangle(mLanNS1( mtalltrials(3,:)>tTH(3) & roisLoc ) , mHanNS1( mtalltrials(3,:)>tTH(3) & roisLoc )  ); % Low - High
% 
%             tonoBNsh2SEL{itroi}(itcon,itrep,:)   =  tonoangle(mLbn2( mtalltrials(1,:)>tTH(1) & roisLoc ) ,   mHbn2( mtalltrials(1,:)>tTH(1) & roisLoc )  ); % Low - High
%             tonoANsh2SEL{itroi}(itcon,itrep,:)   =  tonoangle(mLan2( mtalltrials(2,:)>tTH(2) & roisLoc ) ,   mHan2( mtalltrials(2,:)>tTH(2) & roisLoc )  ); % Low - High
%             tonoANnssh2SEL{itroi}(itcon,itrep,:) =  tonoangle(mLanNS2( mtalltrials(3,:)>tTH(3) & roisLoc ) , mHanNS2( mtalltrials(3,:)>tTH(3) & roisLoc )  ); % Low - High
       

            prefBN{itroi,itcon}    = abs( zscore( mLbn( mtalltrials(1,:)>tTH(1) & roisLoc )) -  zscore(mHbn( mtalltrials(1,:)>tTH(1) & roisLoc )))./(abs(zscore(mLbn( mtalltrials(1,:)>tTH(1) & roisLoc ))) + abs(zscore(mHbn( mtalltrials(1,:)>tTH(1) & roisLoc ))) );
            prefAN{itroi,itcon}    = abs( zscore( mLan( mtalltrials(1,:)>tTH(1) & roisLoc )) -  zscore(mHan( mtalltrials(1,:)>tTH(1) & roisLoc )))./(abs(zscore(mLan( mtalltrials(1,:)>tTH(1) & roisLoc ))) + abs(zscore(mHan( mtalltrials(1,:)>tTH(1) & roisLoc ))) );
            prefANnn{itroi,itcon}  = abs( zscore( mLanNN( mtalltrials(1,:)>tTH(1) & roisLoc )) -  zscore(mHanNN( mtalltrials(1,:)>tTH(1) & roisLoc )))./(abs(zscore(mLanNN( mtalltrials(1,:)>tTH(1) & roisLoc ))) + abs(zscore(mHanNN( mtalltrials(1,:)>tTH(1) & roisLoc ))) );

            
            
            tonoBNsh1{itroi}(itcon,itrep,:)   = (zscore(mLbn1(roisLoc) )   -  zscore(mHbn1(roisLoc))); % Low - High
            tonoANsh1{itroi}(itcon,itrep,:)   = (zscore(mLan1(roisLoc))    -  zscore(mHan1(roisLoc))); % Low - High
            tonoANnssh1{itroi}(itcon,itrep,:) = (zscore(mLanNS1(roisLoc))  -  zscore(mHanNS1(roisLoc))); % Low - High

            tonoBNsh2{itroi}(itcon,itrep,:)   = (zscore(mLbn2(roisLoc))    -  zscore(mHbn2(roisLoc))); % Low - High
            tonoANsh2{itroi}(itcon,itrep,:)   = (zscore(mLan2(roisLoc))    -  zscore(mHan2(roisLoc))); % Low - High
            tonoANnssh2{itroi}(itcon,itrep,:) = (zscore(mLanNS2(roisLoc))  -  zscore(mHanNS2(roisLoc))); % Low - High            
            
            tonoBNsh1SEL{itroi}(itcon,itrep,:)   = ( zscore(mLbn1( mtalltrials(1,:)>tTH(1) & roisLoc )) - zscore(mHbn1( mtalltrials(1,:)>tTH(1) & roisLoc ))  );%./(abs(zscore(mLbn1( mtalltrials(1,:)>tTH(1) & roisLoc )))        +abs(zscore(mHbn1( mtalltrials(1,:)>tTH(1) & roisLoc   )) )); % Low - High
            tonoANsh1SEL{itroi}(itcon,itrep,:)   = ( zscore(mLan1( mtalltrials(2,:)>tTH(2) & roisLoc )) - zscore(mHan1( mtalltrials(2,:)>tTH(2) & roisLoc ))  );%./(abs(zscore(mLan1( mtalltrials(2,:)>tTH(2) & roisLoc )))        + abs(zscore(mHan1( mtalltrials(2,:)>tTH(2) & roisLoc  )) )); % Low - High
            tonoANnssh1SEL{itroi}(itcon,itrep,:) = ( zscore(mLanNS1( mtalltrials(3,:)>tTH(3) & roisLoc )) - zscore(mHanNS1( mtalltrials(3,:)>tTH(3) & roisLoc ))  );%./(abs(zscore(mLanNS1( mtalltrials(3,:)>tTH(3) & roisLoc ))) + abs(zscore(mHanNS1( mtalltrials(3,:)>tTH(3) & roisLoc )) )); % Low - High

            tonoBNsh2SEL{itroi}(itcon,itrep,:)   = ( zscore(mLbn2( mtalltrials(1,:)>tTH(1) & roisLoc ))   - zscore(mHbn2( mtalltrials(1,:)>tTH(1) & roisLoc ))  );%./sum(abs(zscore(mLbn2( mtalltrials(1,:)>tTH(1) & roisLoc )))+abs(zscore(mHbn2( mtalltrials(1,:)>tTH(1) & roisLoc )) )); % Low - High
            tonoANsh2SEL{itroi}(itcon,itrep,:)   = ( zscore(mLan2( mtalltrials(2,:)>tTH(2) & roisLoc ))   - zscore(mHan2( mtalltrials(2,:)>tTH(2) & roisLoc ))  );%./sum(abs(zscore(mLan2( mtalltrials(2,:)>tTH(2) & roisLoc )))+abs(zscore(mHan2( mtalltrials(2,:)>tTH(2) & roisLoc )))); % Low - High
            tonoANnssh2SEL{itroi}(itcon,itrep,:) = ( zscore(mLanNS2( mtalltrials(3,:)>tTH(3) & roisLoc )) - zscore(mHanNS2( mtalltrials(3,:)>tTH(3) & roisLoc ))  );%./sum(abs(zscore(mLanNS2( mtalltrials(3,:)>tTH(3) & roisLoc )))+abs(zscore(mHanNS2( mtalltrials(3,:)>tTH(3) & roisLoc )))); % Low - High


%             tonoBNsh1SEL{itroi}(itcon,itrep,:)   = abs( zscore(mLbn1( mtalltrials(1,:)>tTH(1) & roisLoc )) - zscore(mHbn1( mtalltrials(1,:)>tTH(1) & roisLoc ))  )./(abs(zscore(mLbn1( mtalltrials(1,:)>tTH(1) & roisLoc )))        +abs(zscore(mHbn1( mtalltrials(1,:)>tTH(1) & roisLoc   )) )); % Low - High
%             tonoANsh1SEL{itroi}(itcon,itrep,:)   = abs( zscore(mLan1( mtalltrials(2,:)>tTH(2) & roisLoc )) - zscore(mHan1( mtalltrials(2,:)>tTH(2) & roisLoc ))  )./(abs(zscore(mLan1( mtalltrials(2,:)>tTH(2) & roisLoc )))        + abs(zscore(mHan1( mtalltrials(2,:)>tTH(2) & roisLoc  )) )); % Low - High
%             tonoANnssh1SEL{itroi}(itcon,itrep,:) = abs( zscore(mLanNS1( mtalltrials(3,:)>tTH(3) & roisLoc )) - zscore(mHanNS1( mtalltrials(3,:)>tTH(3) & roisLoc ))  )./(abs(zscore(mLanNS1( mtalltrials(3,:)>tTH(3) & roisLoc ))) + abs(zscore(mHanNS1( mtalltrials(3,:)>tTH(3) & roisLoc )) )); % Low - High
% 
%             tonoBNsh2SEL{itroi}(itcon,itrep,:)   = abs( zscore(mLbn2( mtalltrials(1,:)>tTH(1) & roisLoc ))   - zscore(mHbn2( mtalltrials(1,:)>tTH(1) & roisLoc ))  )./sum(abs(zscore(mLbn2( mtalltrials(1,:)>tTH(1) & roisLoc )))+abs(zscore(mHbn2( mtalltrials(1,:)>tTH(1) & roisLoc )) )); % Low - High
%             tonoANsh2SEL{itroi}(itcon,itrep,:)   = abs( zscore(mLan2( mtalltrials(2,:)>tTH(2) & roisLoc ))   - zscore(mHan2( mtalltrials(2,:)>tTH(2) & roisLoc ))  )./sum(abs(zscore(mLan2( mtalltrials(2,:)>tTH(2) & roisLoc )))+abs(zscore(mHan2( mtalltrials(2,:)>tTH(2) & roisLoc )))); % Low - High
%             tonoANnssh2SEL{itroi}(itcon,itrep,:) = abs( zscore(mLanNS2( mtalltrials(3,:)>tTH(3) & roisLoc )) - zscore(mHanNS2( mtalltrials(3,:)>tTH(3) & roisLoc ))  )./sum(abs(zscore(mLanNS2( mtalltrials(3,:)>tTH(3) & roisLoc )))+abs(zscore(mHanNS2( mtalltrials(3,:)>tTH(3) & roisLoc )))); % Low - High


        end
    end
end
%% RSA matrices Cond x Cond: of split half correlation --> this is figure that we leave out 
% RSA matrices between conditions averaged across split half repetitions
% for each ROI
%  rsaBN   = double(get_RSA_tono(tonoBN));
%  rsaAN   = double(get_RSA_tono(tonoAN));
%  rsaANnn = double(get_RSA_tono(tonoANnn));
%  
%  
  rsaBN   = get_RSA_tono_whithin(tonoBNsh1SEL,tonoBNsh2SEL);
  rsaAN   = get_RSA_tono_whithin(tonoANsh1SEL,tonoANsh2SEL);
  rsaANnn = get_RSA_tono_whithin(tonoANnssh1SEL,tonoANnssh2SEL);
  figure('Color',[1,1,1],'Position',[281 19 1160 820]);
  tonocon = {'Comp','OddB','Unexp'};
  subplot(231);imagesc(squeeze(rsaBN(3,:,:)));  colorbar;title(['Before: Selectivity-RSA (split half) '],'interpreter','Latex'); cl = caxis();
  h = gca(); set(h,'XTick',1:3);set(h,'YTick',1:3);set(h,'XTickLabels',tonocon);set(h,'YTickLabels',tonocon);
  subplot(232);imagesc(squeeze(rsaAN(3,:,:)));  colorbar; title(['After: Selectivity RSA (split half) '],'interpreter','Latex');       caxis(cl);
  h = gca(); set(h,'XTick',1:3);set(h,'YTick',1:3);set(h,'XTickLabels',tonocon);set(h,'YTickLabels',tonocon);
  subplot(233);imagesc(squeeze(rsaANnn(3,:,:)));colorbar;title(['AfterNN: Selectivity RSA (split half) '],'interpreter','Latex');     caxis(cl);
  h = gca(); set(h,'XTick',1:3);set(h,'YTick',1:3);set(h,'XTickLabels',tonocon);set(h,'YTickLabels',tonocon);

  
  rsaBNAN     = double(get_RSA_tono(tonoBN,tonoAN));
  rsaBNANnn   = double(get_RSA_tono(tonoBN,tonoANnn));
  rsaANANnn   = double(get_RSA_tono(tonoAN,tonoANnn));


%%
% figure('Color',[1,1,1]);
figure('Color',[1,1,1],'Position',[281 19 1160 820]),
for itroi = 1:nroi
    
    subplot(nroi, 2, 2*(itroi-1)+1);
    tonoOverlap    = [sum(sign(tonoBNsh1{itroi}) ==  sign(tonoBNsh2{itroi}),3)/nvoxROI(itroi)  ;  sum(sign(tonoANsh1{itroi}) ==  sign(tonoANsh2{itroi}),3)/nvoxROI(itroi) ];
    tonoOverlap    = squeeze(tonoOverlap);
%     boxplot(tonoOverlap')
%     bar([mean(tonoOverlap(1:3,:),2)  mean(tonoOverlap(4:6,2),2)])
    methID = repmat([1 1 1 2 2 2]',1,nrep);
    condID = repmat([1 2 3 1 2 3 ]',1,nrep);
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
    tonoOverlapSel = [ sum(sign(tonoBNsh1SEL{itroi}) ==  sign(tonoBNsh2SEL{itroi}), 3) /size(tonoBNsh1SEL{itroi},3)  ;  sum(sign(tonoANsh1SEL{itroi}) ==  sign(tonoANsh2SEL{itroi}), 3) /size(tonoANsh1SEL{itroi},3)   ];
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

h = gcf;
saveas(h,[folderOut,'TonotopyReliability_withANDwithoutSel_perROI.fig']);
save([folderOut,'TonotopyReliability_withANDwithoutSel_perROI.mat'],'tonoOverlap','tonoBNsh1','tonoBNsh2', 'tonoANsh1', 'tonoANsh2', 'methID', 'condID','nvoxROI', 'ROIname', 'tonoOverlapSel','tonoBNsh1SEL', 'tonoBNsh2SEL', 'tonoANsh1SEL', 'tonoANsh2SEL');
%% Ignore Oddball and Unexp,focus on Complete
condsel  = [1 4 7]; %Complete is in position 1 and 4 of the tonoOverlap matrix
condsel1 = 1;
roi      = zeros(nroi ,3, nrep);
alpha    = 2;
clear tonoOverlap; clear tonoOverlapSel;
for itroi = 1:nroi
%     tonoOverlap    = [sum(tonoBNsh1{itroi} ==  tonoBNsh2{itroi},3)/nvoxROI(itroi)  ;  sum(tonoANsh1{itroi} ==  tonoANsh2{itroi},3)/nvoxROI(itroi); sum(tonoANnssh1{itroi} ==  tonoANnssh2{itroi},3)/nvoxROI(itroi) ];
     tonoOverlap1    = [sum(tonoBNsh1{itroi} ==  tonoBNsh2{itroi},3)/nvoxROI(itroi)  ;  sum(tonoANsh1{itroi} ==  tonoANsh2{itroi},3)/nvoxROI(itroi); sum(tonoANnssh1{itroi} ==  tonoANnssh2{itroi},3)/nvoxROI(itroi) ];
     tonoOverlapROI1(itroi,:,:) = tonoOverlap1(condsel,:);

    tonoOverlap(1,:)    =  columnwisecorrelation(squeeze(tonoBNsh1{itroi}(condsel1,:,:))'  , squeeze(tonoBNsh2{itroi}(condsel1,:,:))','pearson' );
    tonoOverlap(2,:)    =  columnwisecorrelation(squeeze(tonoANsh1{itroi}(condsel1,:,:))'  , squeeze(tonoANsh2{itroi}(condsel1,:,:))','pearson' );
    tonoOverlap(3,:)    =  columnwisecorrelation(squeeze(tonoANnssh1{itroi}(condsel1,:,:))'  , squeeze(tonoANnssh2{itroi}(condsel1,:,:))','pearson' );
    tonoOverlapROI(itroi,:,:) = tonoOverlap;

    tonoOverlapSel1 = [ sum(tonoBNsh1SEL{itroi} ==  tonoBNsh2SEL{itroi}, 3) /size(tonoBNsh1SEL{itroi},3)   ; sum(tonoANsh1SEL{itroi} ==  tonoANsh2SEL{itroi}, 3) /size(tonoANsh1SEL{itroi},3) ; sum(tonoANnssh1SEL{itroi} ==  tonoANnssh2SEL{itroi}, 3) /size(tonoANnssh1SEL{itroi},3)   ];
    tonoOverlapSelROI1(itroi,:,:) = tonoOverlapSel1(condsel,:);
    
    tonoOverlapSel(1,:)    =  columnwisecorrelation(squeeze(tonoBNsh1SEL{itroi}(condsel1,:,:))'  , squeeze(tonoBNsh2SEL{itroi}(condsel1,:,:))','pearson' );
    tonoOverlapSel(2,:)    =  columnwisecorrelation(squeeze(tonoANsh1SEL{itroi}(condsel1,:,:))'  , squeeze(tonoANsh2SEL{itroi}(condsel1,:,:))','pearson' );
    tonoOverlapSel(3,:)    =  columnwisecorrelation(squeeze(tonoANnssh1SEL{itroi}(condsel1,:,:))'  , squeeze(tonoANnssh2SEL{itroi}(condsel1,:,:))','pearson' );
    tonoOverlapSelROI(itroi,:,:) = tonoOverlapSel;

    
    
    roi(itroi,:,:)               = alpha.*itroi;
%                      temp       = corr([ squeeze( mean(tonoBNsh1{itroi},2 )) ; squeeze( mean(tonoBNsh2{itroi},2 )) ]' );
%    corrBN(itroi)
end

meth = 0.*tonoOverlapROI; meth(:,1,:)=1; meth(:,2,:)=2;meth(:,3,:)=3;

figure('Color',[1,1,1],'Position',[281 19 1160 820]),
subplot(221);
boxchart(roi(:),tonoOverlapROI(:),'GroupByColor',meth(:))
% ylim([0.45 0.6])
h = gca; set(h,'FontSize',12);set(h,'XTick',squeeze(roi(:,1,1)))
set(h,'XTickLabel',ROIname,'TickLabelInterpreter','Latex')
ylabel('Consistency','Interpreter','Latex');
title(['no selection '],'interpreter','Latex')


subplot(222);
boxchart(roi(:),tonoOverlapSelROI(:),'GroupByColor',meth(:))
% ylim([0.45 0.6])
h = gca; set(h,'FontSize',12);set(h,'XTick',squeeze(roi(:,1,1)))
set(h,'XTickLabel',ROIname,'TickLabelInterpreter','Latex')
ylabel('Consistency','Interpreter','Latex')
title([' selection with F p0 = ', num2str(p0) ],'interpreter','Latex')

h = gcf;
saveas(h,[folderOut,'TonotopyReliability_CompleteOnly.fig']);
save([folderOut,'TonotopyReliability_CompleteOnly.mat'],'tonoOverlap','roi', 'condsel','condsel1', 'roi','alpha', 'tonoOverlap1', 'tonoBNsh1', 'tonoBNsh2', 'tonoANsh1', 'tonoANsh2', 'tonoOverlapROI1', 'tonoOverlapSel1', 'tonoOverlapSelROI1', 'tonoOverlapSel', 'meth','nvoxROI', 'ROIname','tonoBNsh1SEL', 'tonoBNsh2SEL', 'tonoANsh1SEL', 'tonoANsh2SEL');

