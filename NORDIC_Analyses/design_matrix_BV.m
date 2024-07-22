function [X1 ,prednames] = design_matrix_BV(X,nvol,TR)
% creates a design matrix identical of BV
pttp = 5;
nttp = 16;
pnr = 6;
ons = 0;
pdsp = 1;
ndsp = 1;
nTP  = nvol;
prtr = TR;
params = struct('erlen',0,'hshape','twogamma','hpttp',pttp,'hnttp',nttp,'hpnr',pnr,'hons',ons,'hpdsp',pdsp,'hndsp',ndsp,'nvol',nTP,'prtr',prtr,'rcond',X.NrOfConditions+1);
X1 = X.CreateSDM(params); 
prednames = X1.PredictorNames;
X1 = X1.SDMMatrix;
