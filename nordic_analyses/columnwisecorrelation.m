function [corrAB ] = columnwisecorrelation(A,B,type)
%% column by column correlation between two matrices
switch type
    case 'pearson'
%% column by column correlation between two matrices
An=bsxfun(@minus,A,mean(A,1)); %%% zero-mean
Bn=bsxfun(@minus,B,mean(B,1)); %%% zero-mean
% covAB = sum(An.*Bn,1)./(size(A,1)-1);
An=bsxfun(@times,An,1./sqrt(sum(An.^2,1))); %% L2-normalization
Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1))); %% L2-normalization
corrAB =sum(An.*Bn,1);
    case 'spearman'
        corrAB = 0.*A(1,:);
        for it = 1:size(A,2)
          loc  = ~isnan(A(:,it)) & ~isnan(B(:,it));
          corrAB(it) = corr(A(loc,it),B(loc,it),'type','spearman');  
        end
end