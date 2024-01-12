function C = get_RSA_tono(tono1,tono2)
% RSA matrices between conditions averaged across split half repetitions
% for each ROI
% for itROI = 1:size(tonoSH1,2)
%       nrep  = size(tonoSH1{itROI},2);
%       for itrep = 1:nrep
%           C(itROI,itrep,:,:) = corr( squeeze(tonoSH1{itROI}(:,itrep,:))' ,  squeeze(tonoSH2{itROI}(:,itrep,:))' );
% %          C(itROI,itrep,:,:) = corr( squeeze(tonoSH1{itROI}(:,itrep,:))' ,  squeeze(tonoSH2{itROI}(:,itrep,:))' );
% % 
%      end
% %             t1      = squeeze(mean( squeeze(tonoSH1{itROI}),2 ));
% %             t2      = squeeze(mean( squeeze(tonoSH2{itROI}),2 ));
% %      for itcond1 = 1:size(t1,1)
% %          for itcond2 = 1:size(t2,1)
% % 
% %           C(itROI,itcond1,itcond2) = corr(t1(itcond1,:)' , t2(itcond2,:)'  );
% % 
% % %           C(itROI,itcond1,itcond2) = corr( squeeze(tonoSH1{itcond1,itROI}') , squeeze(tonoSH2{itcond2,itROI}' ) );
% % %           C(itROI,itcond1,itcond2) = sum( squeeze(tonoSH1{itcond1,itROI}') == squeeze(tonoSH2{itcond2,itROI}' ) )/length( squeeze(tonoSH2{itcond2,itROI}' ) );        
% %          end
% %      end
% end
% C = squeeze(mean(C,2));
% disp('A')
for itcond1 = 1:size(tono1,1)
    for itcond2 = 1:size(tono1,1)
        for itroi = 1:size(tono1,2)
            C(itroi,itcond1 ,itcond2 ) = corr( tono1{itcond1,itroi}'  ,tono2{itcond2,itroi}' );
%             C(itcond1 ,itcond2, itroi) = sum( tono{itcond1,itroi}'  == tono{itcond2,itroi}' )/length(tono{itcond1,itroi});
        end
    end
end