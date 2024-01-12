function [Fbn,Fan,FanNS] = SaveFMaps(connames,colnames,Xbdiag,bBn,bAn,bAnNS,tsbn,tsan,tsanNS,ccInv,folderOut,roi_large,VTCInfo, subj)

colnames = string(colnames);
connames = string(connames);
contrast = contains(colnames, connames); %%% for F Map - All 1's except for the constant term.
tsbn_hat = Xbdiag*bBn; e_bn = tsbn - tsbn_hat; ve_bn = var(e_bn); Fbn = (bBn'*contrast')./(sqrt(ve_bn'*contrast*ccInv*contrast')); % For Fbn see formula t in notes
tsan_hat = Xbdiag*bAn; e_an = tsan - tsan_hat; ve_an = var(e_an); Fan = (bAn'*contrast')./(sqrt(ve_an'*contrast*ccInv*contrast'));
tsanNS_hat = Xbdiag*bAnNS; e_anNS = tsanNS - tsanNS_hat; ve_anNS = var(e_anNS); FanNS = (bAnNS'*contrast')./(sqrt(ve_anNS'*contrast*ccInv*contrast'));
create_vmp(folderOut,[subj,'_FMaps.vmp'],VTCInfo,[Fbn,Fan,FanNS],roi_large,{'Fbn','Fan','Fan_NoNoise'});

A = bBn(contrast,:);
tst_bn =sqrt(size(A,1)).* (mean(A)./std(A));
p_bn = 2*tcdf(abs(tst_bn),size(A,1)-1,'upper');
[fdr,q] = mafdr(p_bn);
index_y = find(q<0.01);
tst_bn_thr = zeros(size(tst_bn));
tst_bn_thr(index_y) = tst_bn(index_y);
disp(['Before Nordic ' , num2str(length(index_y)/length(tst_bn) ) , ' significant voxels p90_t ' , num2str(prctile(tst_bn(index_y),90) ) ])

A = bAn(contrast,:);
tst_An = sqrt(size(A,1)).* (mean(A)./std(A));
p_An = 2*tcdf(abs(tst_An),size(A,1)-1,'upper');
[fdr,q] = mafdr(p_An);
index_y = find(q<0.01);
tst_An_thr = zeros(size(tst_An));
tst_An_thr(index_y) = tst_An(index_y);
disp(['Nordic ' , num2str(length(index_y)/length(tst_An)) , ' significant voxels p90_t ' , num2str(prctile(tst_An(index_y),90) ) ])


A = bAnNS(contrast,:);
tst_AnNS = sqrt(size(A,1)).* (mean(A)./std(A));
p_AnNS = 2*tcdf(abs(tst_AnNS),size(A,1)-1,'upper');
[fdr,q] = mafdr(p_AnNS);
index_y = find(q<0.01);
tst_AnNS_thr = zeros(size(tst_AnNS));
tst_AnNS_thr(index_y) = tst_AnNS(index_y);
disp(['Nordic NN ' , num2str(length(index_y)/length(tst_AnNS)) , ' significant voxels p90_t ' , num2str(prctile(tst_AnNS(index_y),90) ) ])


create_vmp(folderOut, [subj, '_FMaps_singleTrials.vmp'],VTCInfo,[tst_bn',tst_An',tst_AnNS', tst_bn_thr', tst_An_thr', tst_AnNS_thr'],roi_large,{'Fbn','Fan','Fan_NoNoise','Fbn_th','Fan_th','Fan_NoNoise_th'});


%% For tonotopy

A = bBn(contrast,:);
tst_bn =sqrt(size(A,1)).* (mean(A)./std(A));
p_bn = 2*tcdf(abs(tst_bn),size(A,1)-1,'upper');
[fdr,q] = mafdr(p_bn);
index_y = find(q<0.0001);
tst_bn_thr = zeros(size(tst_bn));
tst_bn_thr(index_y) = tst_bn(index_y);


A = bAn(contrast,:);
tst_An = sqrt(size(A,1)).* (mean(A)./std(A));
p_An = 2*tcdf(abs(tst_An),size(A,1)-1,'upper');
[fdr,q] = mafdr(p_An);
index_y = find(q<0.0001);
tst_An_thr = zeros(size(tst_An));
tst_An_thr(index_y) = tst_An(index_y);


A = bAnNS(contrast,:);
tst_AnNS = sqrt(size(A,1)).* (mean(A)./std(A));
p_AnNS = 2*tcdf(abs(tst_AnNS),size(A,1)-1,'upper');
[fdr,q] = mafdr(p_AnNS);
index_y = find(q<0.0001);
tst_AnNS_thr = zeros(size(tst_AnNS));
tst_AnNS_thr(index_y) = tst_AnNS(index_y);


create_vmp(folderOut, [subj, '_FMaps_singleTrials_ForTonotopy.vmp'],VTCInfo,[tst_bn_thr', tst_An_thr', tst_AnNS_thr'],roi_large,{'Fbn_th','Fan_th','Fan_NoNoise_th'});