function [TSS,SSb,SSe] = partitionSSols(Y,X,beta)
TSS = sum( (Y - mean(Y) ).^2  ); % total sum of squares
e   = Y - X*beta;                % residuals
SSb = sum( (X*beta - mean(X*beta)).^2  );
SSe = sum( (e - mean(e) ).^2  );
