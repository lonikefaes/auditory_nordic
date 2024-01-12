function [CorrBetas, tonoBn, tonoAn, tonoAnNS] = CorrBetasAcrossCond(connames,colnames,bBn,bAn,bAnNS,th,folderOut,roi_large,VTCInfo,subj, excludevox)

Bn = [];
An = [];
AnNS = [];

Bn2 = {};
An2 = {};
AnNS2 = {};
for i=1:size(connames,1)
    cons = contains(colnames, {'Constant'});
    cons = find(cons==1);

    for j=1:size(connames,2)
        select = contains(colnames, connames{i,j});
        index = find(select==1);
        for z = 1:length(cons)
            if z==1
                indexr = find(index<cons(z));
            else
                indexr = find(index>cons(z-1) & index<cons(z));
            end

            Bn(i,j,z,:)   = mean(bBn(index(indexr),:));
            An(i,j,z,:)   = mean(bAn(index(indexr),:));
            AnNS(i,j,z,:) = mean(bAnNS(index(indexr),:));

            Bn2{i,j,z}   = (bBn(index(indexr),:));
            An2{i,j,z}   = (bAn(index(indexr),:));
            AnNS2{i,j,z} = (bAnNS(index(indexr),:));
        end 
    end
end


Bn   = reshape(Bn, size(Bn,1)*size(Bn,2),size(Bn,3),size(Bn,4));
An   = reshape(An, size(An,1)*size(An,2),size(An,3), size(An,4));
AnNS = reshape(AnNS, size(AnNS,1)*size(AnNS,2),size(AnNS,3), size(AnNS,4));

for j=1:size(Bn,2)
    for i=1:size(Bn,1)
        CorrBetas(1,i,j) = corr2(squeeze(Bn(i,j,:)),squeeze(An(i,j,:)));
        CorrBetas(2,i,j) = corr2(squeeze(Bn(i,j,:)),squeeze(AnNS(i,j,:)));
    
    end
end 


%% Tonotopy %%%
%%% take high and low complete conditions
Bn   = Bn([1,4],:,:);
An   = An([1,4],:,:);
AnNS = AnNS([1,4],:,:);
 
%% threshold - we use a threshold computed on all the runs (this way the consistency is not dependent on the variability of selected voxels)
CompH = [];
CompL = [];

for i=1:size(Bn2,3)
    CompH = [CompH;Bn2{1,1,i}];
    CompL = [CompL;Bn2{1,2,i}];
end

All  = [CompH;CompL];
F = sqrt(size(All,1)).*(mean(All)./std(All));

th = 1; % thres is coming from outside 1 means no selection
th = 0.05;
[sortF,indexSort] = sort(F,'descend');
%%% best 80% voxels
nrvox  = round(th*length(sortF));
selvox = indexSort(1:nrvox);

%%
tt = squeeze(mean(Bn,2));
tt2 = squeeze(mean(An,2));
tt3 = squeeze(mean(AnNS,2));
tt_z = zscore(tt,[],2);
tt2_z = zscore(tt2,[],2);
tt3_z = zscore(tt3,[],2);

[mm,tono]  = max(tt_z);
tonoAverageBN = tono;% tonotopy from the average across runs in original data
[mm,tono2] = max(tt2);
[mm,tono3] = max(tt3);

tono(find(tono==1)) = -10;
tono(find(tono==2)) = 10;

tono2(find(tono2==1)) = -10;
tono2(find(tono2==2)) = 10;

tono3(find(tono3==1)) = -10;
tono3(find(tono3==2)) = 10;

Map = [tono]';

create_vmp_tono(folderOut, [subj, '_TonoNoSel.vmp'],VTCInfo,Map,roi_large,{'tono'});


%% Do tonotopy per run and correlate NORDIC (NN) vs original
for i=1:size(Bn,2)
    tempBn   = zscore(squeeze(Bn(:,i,:)),[],2);
    tempAn   = zscore(squeeze(An(:,i,:)),[],2);
    tempAnNS = zscore(squeeze(AnNS(:,i,:)),[],2);

%     tempBn   = (squeeze(Bn(:,i,:))  );
%     tempAn   = (squeeze(An(:,i,:))  );
%     tempAnNS = (squeeze(AnNS(:,i,:)));

    [mm, t] = max(tempBn);
    tonoBn(:,i) = t;
    [mm, t] = max(tempAn);
    tonoAn(:,i) = t;
    [mm,t] = max(tempAnNS);
    tonoAnNS(:,i) = t;

%     tonoBn(:,i)   = tempBn;
%     tonoAn(:,i)   = tempAn;
%     tonoAnNS(:,i) = tempAnNS;


    corrTono2Orig_N(i)  =  corr(tonoBn(selvox,i),tonoAn(selvox,i));
    corrTono2Orig_NN(i) =  corr(tonoBn(selvox,i),tonoAnNS(selvox,i));

    corrTono2Orig_Nave(i)  =  corr(tonoAverageBN(selvox)', tonoAn(selvox,i));      % Average tono of BN correlated with each run of NORDIC
    corrTono2Orig_NNave(i) =  corr(tonoAverageBN(selvox)', tonoAnNS(selvox,i));

     excludeRun                =  setdiff(1:size(Bn,2) , i );%indices excluding a run
%      excludeRun                = 1:8;
    [~,tonoAVEloroBN]          =  max(zscore(squeeze( mean(Bn(:,excludeRun,:),2) ),[],2) );% average tonotopy with run i left out for consistency check we did: [~,tonoAVE]       = max(zscore(squeeze( mean(Bn,2) ),[],2) ) == tonoAverageBN 
    corrTono2Orig_aveLORO(i)   =  corr(tonoAVEloroBN(selvox)', tonoBn(selvox,i));         % correlation run left out to the average except this run: all within original data
    corrTono2Orig_NaveLORO(i)  =  corr(tonoAVEloroBN(selvox)', tonoAn(selvox,i));         % Average tono of BN correlated with each run of NORDIC
    corrTono2Orig_NNaveLORO(i) =  corr(tonoAVEloroBN(selvox)', tonoAnNS(selvox,i));
end


for i= 1:size(Bn,2)
     excludeRun                 =  setdiff(1:size(Bn,2) , i );%indices excluding a run
     corrTono2Orig_aveLORO(i)   =  corr(mean(tonoBn(selvox,excludeRun),2), tonoBn(selvox,i) );         % correlation run left out to the average except this run: all within original data
     corrTono2Orig_NaveLORO(i)  =  corr(mean(tonoBn(selvox,excludeRun),2), tonoAn(selvox,i) );         % Average tono of BN correlated with each run of NORDIC
     corrTono2Orig_NNaveLORO(i) =  corr(mean(tonoBn(selvox,excludeRun),2), tonoAnNS(selvox,i));

%      corrTono2Orig_aveLORO(i)   =  corr(diff(squeeze(mean(Bn(:,excludeRun,selvox),2)))'  ,     diff(squeeze(Bn(:,i,selvox)))' );         % correlation run left out to the average except this run: all within original data
%      corrTono2Orig_NaveLORO(i)  =  corr(diff(squeeze(mean(Bn(:,excludeRun,selvox),2)))' ,     diff(squeeze(An(:,i,selvox)))'  );         % Average tono of BN correlated with each run of NORDIC
%      corrTono2Orig_NNaveLORO(i) =  corr(diff(squeeze(mean(Bn(:,excludeRun,selvox),2)))' ,     diff(squeeze(AnNS(:,i,selvox)))'  );
end
%      excludeRun1 = [1:4]; excludeRun2 = [5:8];
%      Bn = zscore(Bn,[],2); An = zscore(An,[],2); AnNS = zscore(AnNS,[],2);
%      corrTono2Orig_aveLORO(i)   =  corr(diff(squeeze(mean(Bn(:,excludeRun1,selvox),2)))'  ,    diff(squeeze(mean(Bn(:,excludeRun2,selvox),2)))'  );         % correlation run left out to the average except this run: all within original data
%      corrTono2Orig_NaveLORO(i)  =  corr(diff(squeeze(mean(Bn(:,excludeRun1,selvox),2)))' ,     diff(squeeze(mean(An(:,excludeRun2,selvox),2)))'   );         % Average tono of BN correlated with each run of NORDIC
%      corrTono2Orig_NNaveLORO(i) =  corr(diff(squeeze(mean(Bn(:,excludeRun1,selvox),2)))' ,     diff(squeeze(mean(AnNS(:,excludeRun2,selvox),2)))'   );


% methCount = [repmat(1,size(Bn,2),1); repmat(2,size(Bn,2),1)];
methCount = [repmat(1,size(Bn,2),1); repmat(2,size(Bn,2),1); repmat(3,size(Bn,2),1) ];



%% equivalent figure with the average across runs for the BN tonotopy
figure
% mb = boxchart([corrTono2Orig_Nave';corrTono2Orig_NNave'], 'GroupByColor', methCount, 'WhiskerLineStyle' ,'none','MarkerStyle','none');
% mb = boxchart([corrTono2Orig_aveLORO'; corrTono2Orig_Nave';corrTono2Orig_NNave'], 'GroupByColor', methCount, 'WhiskerLineStyle' ,'none','MarkerStyle','none');
mb = boxchart([corrTono2Orig_aveLORO'; corrTono2Orig_NaveLORO';corrTono2Orig_NNaveLORO'], 'GroupByColor', methCount, 'WhiskerLineStyle' ,'none','MarkerStyle','none');

mb(2).BoxFaceColor = [144 103 167]./256;
mb(3).BoxFaceColor = [0.4660 0.6740 0.1880];
title('Run to run correlation of tonotopic maps to average tonotopy BN')
xlabel('Method')
ylabel('Correlation')
set(gca, 'XTick', [])
% legend('NORDIC vs original', 'NORDIC NN vs original')
  legend( 'loro ORIG' ,'NORDIC vs original', 'NORDIC NN vs original')


s = gcf;

%saveas(s,[folderOut, [subj,'_CorrelationTonotopy_run2Average.fig']]);
%save([folderOut,[subj, '_CorrelationTonotopy_run2Average.mat']],'corrTono2Orig_Nave', 'corrTono2Orig_NNave');
%saveas(s, [folderOut, [subj,'_CorrelationTonotopy_run2Average.svg']]);


%% Pairwise correlations of the betas

% Corr at the beta level between run
PPCor_Bn   = corr(  squeeze(Bn(1,:,selvox)  )'  );
PPCor_An   = corr(  squeeze(An(1,:,selvox)  )'  );
PPCor_AnNS = corr(  squeeze(AnNS(1,:,selvox))'  );

PPCor_Bnave   = corr(  squeeze(Bn(1,:,selvox))' ,   squeeze(mean(Bn(1,:,selvox))   )  );
PPCor_Anave   = corr(  squeeze(An(1,:,selvox))' ,   squeeze(mean(An(1,:,selvox))   )  );
PPCor_AnNSave = corr(  squeeze(AnNS(1,:,selvox))' , squeeze(mean(AnNS(1,:,selvox))   )  );


PPCor_BnaveBN   = corr(  squeeze(Bn(1,:,selvox))' ,   squeeze(mean(Bn(1,:,selvox))   )  );
PPCor_AnaveBN   = corr(  squeeze(An(1,:,selvox))' ,   squeeze(mean(Bn(1,:,selvox))   )  );
PPCor_AnNSaveBN = corr(  squeeze(AnNS(1,:,selvox))' , squeeze(mean(Bn(1,:,selvox))   )  );

for i = 1:size(Bn,2)

    excludeRun         =  setdiff(1:size(Bn,2) , i );%indices excluding a run
    PPCor_BnaveBN(i)   = corr(  squeeze(Bn(1,i,selvox)) ,   squeeze(mean(Bn(1,excludeRun,selvox))   )  );
    PPCor_AnaveBN(i)   = corr(  squeeze(An(1,i,selvox)) ,   squeeze(mean(Bn(1,excludeRun,selvox))   )  );
    PPCor_AnNSaveBN(i) = corr(  squeeze(AnNS(1,i,selvox)) , squeeze(mean(Bn(1,excludeRun,selvox))   )  );
end


PPCor_Bn = triu(PPCor_Bn,1);
PPCor_An = triu(PPCor_An,1);
PPCor_AnNS = triu(PPCor_AnNS,1);
PPCor_Bn = PPCor_Bn(find(PPCor_Bn~=0));
PPCor_An = PPCor_An(find(PPCor_An~=0));
PPCor_AnNS = PPCor_AnNS(find(PPCor_AnNS~=0));

methCount = [repmat(1,size(PPCor_Bn,1),1); repmat(2,size(PPCor_Bn, 1),1); repmat(3,size(PPCor_Bn, 1),1)];

figure
boxchart([PPCor_Bn; PPCor_An; PPCor_AnNS], 'GroupByColor', methCount, 'WhiskerLineStyle' ,'none','MarkerStyle','none')
title('Run to run correlation for each method')
xlabel('Method')
ylabel('Correlation')
set(gca, 'XTick', [])
legend('Original','NORDIC', 'NORDIC NN')


s = gcf;

saveas(s,[folderOut, [subj,'_Run2runCorrelation_perMethod.fig']]);
save([folderOut,[subj, '_Run2runCorrelation_perMethod.mat']],'PPCor_Bn', 'PPCor_An', 'PPCor_AnNS');
saveas(s, [folderOut, [subj,'_Run2runCorrelation_perMethod.svg']]);


%% equivalent figure by comparing betas with the average BETAS of BN

% methCount = [repmat(1,size(PPCor_Bnave, 1),1); repmat(2,size(PPCor_Bnave, 1),1)];
  methCount = [repmat(1,size(PPCor_Bnave, 1),1); repmat(2,size(PPCor_Bnave, 1),1); repmat(3,size(PPCor_Bnave, 1),1)];

figure
% mb = boxchart([   PPCor_AnaveBN; PPCor_AnNSaveBN], 'GroupByColor', methCount, 'WhiskerLineStyle' ,'none','MarkerStyle','none');
mb = boxchart([ PPCor_BnaveBN;  PPCor_AnaveBN; PPCor_AnNSaveBN], 'GroupByColor', methCount, 'WhiskerLineStyle' ,'none','MarkerStyle','none');

mb(2).BoxFaceColor = [144 103 167]./256;
mb(3).BoxFaceColor = [0.4660 0.6740 0.1880];
title('Correlation of betas of each run to average of BN')
xlabel('Method')
ylabel('Correlation')
set(gca, 'XTick', [])
legend('NORDIC vs original', 'NORDIC NN vs original')

s = gcf;

saveas(s,[folderOut, [subj,'_Run2runCorr_AverageBN.fig']]);
save([folderOut,[subj, '_Run2runCorr_AverageBN.mat']],'PPCor_AnaveBN', 'PPCor_AnNSaveBN');
saveas(s, [folderOut, [subj,'_Run2runCorr_AverageBN.svg']]);
