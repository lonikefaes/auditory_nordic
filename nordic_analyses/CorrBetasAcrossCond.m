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


Bn   = reshape(Bn, size(Bn,1)*size(Bn,2),size(Bn,3),size(Bn,4));    % This matrix contains values of all conditions by 8 runs by all voxels
An   = reshape(An, size(An,1)*size(An,2),size(An,3), size(An,4));
AnNS = reshape(AnNS, size(AnNS,1)*size(AnNS,2),size(AnNS,3), size(AnNS,4));

for j=1:size(Bn,2)
    for i=1:size(Bn,1)
        CorrBetas(1,i,j) = corr2(squeeze(Bn(i,j,:)),squeeze(An(i,j,:)));    % Correlate betas of before NORDIC with After NORDIC
        CorrBetas(2,i,j) = corr2(squeeze(Bn(i,j,:)),squeeze(AnNS(i,j,:)));  % Bn and AnNS
        % Conditions are ordere how they are put in connames, so first 3
        % are H conditions and then second 3 are low conditions
    end
end 


%% Tonotopy %%%
%%% take high and low complete conditions
Bn   = Bn([1,4],:,:);       % We take only conditions 1 and 4 --> CompH and CompL
An   = An([1,4],:,:);
AnNS = AnNS([1,4],:,:);
 
%% threshold - we use a threshold computed on all the runs (this way the consistency is not dependent on the variability of selected voxels)
% currently we don't use a threshold, we only exclude a couple voxels that
% have extremely low tSNR outside of this function.
CompH = [];
CompL = [];

for i=1:size(Bn2,3)
    CompH = [CompH;Bn2{1,1,i}];
    CompL = [CompL;Bn2{1,2,i}];
end

All  = [CompH;CompL];
F = sqrt(size(All,1)).*(mean(All)./std(All));

th = 1; % thres is coming from outside 1 means no selection
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

[mm,tono]  = max(tt_z); % Find the max of the z-scores to see whether response is highest for CompH or CompL
tonoAverageBN = tono;% tonotopy from the average across runs in original data
[mm,tono2] = max(tt2);
[mm,tono3] = max(tt3);

tono(find(tono==1)) = -10; % Create maps with CompH == -10 and CompL == 10
tono(find(tono==2)) = 10;

tono2(find(tono2==1)) = -10;
tono2(find(tono2==2)) = 10;

tono3(find(tono3==1)) = -10;
tono3(find(tono3==2)) = 10;

Map = [tono]';

create_vmp_tono(folderOut, [subj, '_TonoNoSel.vmp'],VTCInfo,Map,roi_large,{'tono'});


%% Do tonotopy per run and correlate NORDIC (NN) vs original
for i=1:size(Bn,2)
    tempBn   = zscore(squeeze(Bn(:,i,:)),[],2); % z-score the betas per run
    tempAn   = zscore(squeeze(An(:,i,:)),[],2);
    tempAnNS = zscore(squeeze(AnNS(:,i,:)),[],2);

    [mm, t] = max(tempBn);  % take max --> preference for CompH or CompL
    tonoBn(:,i) = t;
    [mm, t] = max(tempAn);
    tonoAn(:,i) = t;
    [mm,t] = max(tempAnNS);
    tonoAnNS(:,i) = t;


    corrTono2Orig_N(i)  =  corr(tonoBn(selvox,i),tonoAn(selvox,i)); % Correlation run 1 or Bn to run 1 of AN etc.
    corrTono2Orig_NN(i) =  corr(tonoBn(selvox,i),tonoAnNS(selvox,i));
end

methCount = [repmat(1,size(Bn,2),1); repmat(2,size(Bn,2),1)];


%% Pairwise correlations of the betas

% Correlation at the beta level between runs
PPCor_BnH   = corr(  squeeze(Bn(1,:,selvox)  )'  ); %Matrix of pairwise correlations for Bn CompH
PPCor_AnH   = corr(  squeeze(An(1,:,selvox)  )'  );
PPCor_AnNSH = corr(  squeeze(AnNS(1,:,selvox))'  );

PPCor_BnL   = corr(  squeeze(Bn(2,:,selvox)  )'  ); %Matrix of pairwise correlations for Bn CompL
PPCor_AnL   = corr(  squeeze(An(2,:,selvox)  )'  );
PPCor_AnNSL = corr(  squeeze(AnNS(2,:,selvox))'  );


PPCor_BnH = triu(PPCor_BnH,1);
PPCor_AnH = triu(PPCor_AnH,1);
PPCor_AnNSH = triu(PPCor_AnNSH,1);
PPCor_BnH = PPCor_BnH(find(PPCor_BnH~=0));
PPCor_AnH = PPCor_AnH(find(PPCor_AnH~=0));
PPCor_AnNSH = PPCor_AnNSH(find(PPCor_AnNSH~=0));

methCount = [repmat(1,size(PPCor_BnH,1),1); repmat(2,size(PPCor_BnH, 1),1); repmat(3,size(PPCor_BnH, 1),1)];

figure
boxchart([PPCor_BnH; PPCor_AnH; PPCor_AnNSH], 'GroupByColor', methCount, 'WhiskerLineStyle' ,'none','MarkerStyle','none')
title('Run to run correlation for each method')
xlabel('Method', 'FontName', 'Arial', 'FontSize', 12)
ylabel('Correlation', 'FontName', 'Arial', 'FontSize', 12)
set(gca, 'XTick', [])
legend({'Standard','NORDIC', 'NORDIC NN'},'FontName', 'Arial', 'FontSize', 12)
p = gcf;
% here we compared the three possible pairs by separate t-test 
allData = [PPCor_BnH; PPCor_AnH; PPCor_AnNSH];
[H,P_std_vs_N,CI,STATS_std_vs_N]  = ttest2(allData(methCount==1),  allData(methCount==2)  );     t_std_vs_N   = STATS_std_vs_N.tstat; 
[H,P_std_vs_Nnn,CI,STATS_std_vs_Nnn]  = ttest2(allData(methCount==1),  allData(methCount==3)  ); t_std_vs_Nnn = STATS_std_vs_Nnn.tstat; 
[H,P_N_vs_Nnn,CI,STATS_N_vs_Nnn]  = ttest2(allData(methCount==2),  allData(methCount==3)  );     t_N_vs_Nnn    = STATS_N_vs_Nnn.tstat; 



saveas(p,[folderOut, [subj,'_Run2runCorrelation_perMethod.fig']]);
save([folderOut,[subj, '_Run2runCorrelation_perMethod.mat']],'PPCor_BnH', 'PPCor_AnH', 'PPCor_AnNSH','P_std_vs_N','P_std_vs_Nnn','P_N_vs_Nnn','t_std_vs_N','t_std_vs_Nnn','t_N_vs_Nnn');
saveas(p, [folderOut, [subj,'_Run2runCorrelation_perMethod.svg']]);


%% equivalent figure by comparing betas with the average BETAS of BN - LORO

for i = 1:size(Bn,2)

    excludeRun         =  setdiff(1:size(Bn,2) , i );%indices excluding a run
    PPCor_BnaveBNH(i)   = corr(  squeeze(Bn(1,i,selvox)) ,   squeeze(mean(Bn(1,excludeRun,selvox))   )  );
    PPCor_AnaveBNH(i)   = corr(  squeeze(An(1,i,selvox)) ,   squeeze(mean(Bn(1,excludeRun,selvox))   )  );
    PPCor_AnNSaveBNH(i) = corr(  squeeze(AnNS(1,i,selvox)) , squeeze(mean(Bn(1,excludeRun,selvox))   )  );

    PPCor_BnaveBNL(i)   = corr(  squeeze(Bn(2,i,selvox)) ,   squeeze(mean(Bn(2,excludeRun,selvox))   )  );
    PPCor_AnaveBNL(i)   = corr(  squeeze(An(2,i,selvox)) ,   squeeze(mean(Bn(2,excludeRun,selvox))   )  );
    PPCor_AnNSaveBNL(i) = corr(  squeeze(AnNS(2,i,selvox)) , squeeze(mean(Bn(2,excludeRun,selvox))   )  );
end

CorrH = [PPCor_BnaveBNH;PPCor_AnaveBNH; PPCor_AnNSaveBNH];
CorrL = [PPCor_BnaveBNL;PPCor_AnaveBNL; PPCor_AnNSaveBNL];

temp = zeros(6,z);

i = 1;
for j = 1:3
    temp(i,:) = [CorrH(j,:)];
    i = i+1;
    temp(i,:) = [CorrL(j,:)];
    i = i+1;
end

CorrHL = temp;

% CorrHL = cat(3,CorrH, CorrL);
% CorrHL = permute(CorrHL, [3,1,2]);

methodID = repmat([1 2 3 4 5 6]', [1,z]);
meth = repmat([1 1 2 2 3 3]', [1,z]);

% CorrBetasHL = 0.*CorrHL;  CorrBetasHL(1,:,:)  = 1; CorrBetasHL(2,:,:) = 2;
% CorrBetasHL2 = 0.*CorrHL; CorrBetasHL2(:,1,:) = 1; CorrBetasHL2(:,2,:) = 2; CorrBetasHL2(:,3,:) = 3; 


figure('Color',[1,1,1],'Position',[281 19 1160 820]),

mb = boxchart(methodID(:), CorrHL(:),  'GroupByColor', meth(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');

ylabel('Cross Validated Correlation', 'FontName', 'Arial', 'FontSize', 12)
xticklabels(['PredH'; 'PredL'; 'PredH'; 'PredL'; 'PredH'; 'PredL'])
xticks([0.66, 1.66, 3, 4, 5.33, 6.33])
%xlabel('')
%xticks('');
legend({'Standard', 'NORDIC', 'NORDIC NN'}, 'FontName', 'Arial', 'FontSize', 12)
%ylim([0.0 0.4])

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';
% repeated measures anova: the same run is measured for every cond and
% method
% factorCond = 1+ 1. *( methodID == 1 | methodID == 3  | methodID == 5 );% here we create a factor that codes for PredL and PredH 
% t          = table(CorrHL(:,1),CorrHL(:,2),CorrHL(:,3),CorrHL(:,4),CorrHL(:,5), CorrHL(:,6),CorrHL(:,7),CorrHL(:,8) ,meth(:,1) ,  factorCond(:,1), 'VariableNames', {'run1','run2','run3','run4','run5','run6','run7','run8','method','cond'} );
% run_design = table( [1:size(methodID,2)]' , 'VariableNames',{'runID'});
% rm = fitrm(t,'run1-run8 ~ method + cond +method *  cond', 'WithinDesign' , run_design);
% anovatable = ranova(rm);

factorCond = 1 + 1. *( methodID == 1 | methodID == 3  | methodID == 5 );% here we create a factor that codes for PredL and PredH 
% t          = table(CorrHL(:,1),CorrHL(:,2),CorrHL(:,3),CorrHL(:,4),CorrHL(:,5), CorrHL(:,6),CorrHL(:,7),CorrHL(:,8) ,meth(:,1) ,  factorCond(:,1), 'VariableNames', {'run1','run2','run3','run4','run5','run6','run7','run8','method','cond'} );
tCorrHL    = CorrHL';
t          = table( tCorrHL(:,1),tCorrHL(:,2),tCorrHL(:,3), tCorrHL(:,4),tCorrHL(:,5),tCorrHL(:,6), 'VariableNames',{'c1','c2','c3','c4','c5','c6'}     );
whithin    = table( categorical(meth(:,1) ), categorical(factorCond(:,1)) ,'VariableNames', {'method', 'cond'} );
rm = fitrm(t,'c1-c6 ~ 1', 'WithinDesign' , whithin);
% contrast between methods collapsing for conditions 1 is std 2 is Nordic 3 is Nordic no noise
Cstd_vs_N    = [1  1 -1  -1   0  0 ]' ;
CN_vs_Nnn    = [0  0  1   1  -1 -1 ]' ;
Cstd_vs_Nnn  = [1  1  0   0  -1 -1 ]' ;

anovaTable_std_vs_N  = ranova(rm,'WithinModel',Cstd_vs_N);
anovaTable_N_vs_Nss  = ranova(rm,'WithinModel',CN_vs_Nnn);
anovaTable_std_vs_Nnn = ranova(rm,'WithinModel',Cstd_vs_Nnn);
% Contrast between conditions: collapsing for methods
Cpred_vs_unpred            = [1  -1   1   -1   1 -1 ]' ;
anovaTable_pred_vs_unpred  = ranova(rm,'WithinModel',Cpred_vs_unpred);

s = gcf;
saveas(s,[folderOut, [subj,'_SpatialCorrelation_BetasHL.fig']]);
save([folderOut,[subj, '_SpatialCorrelation_BetasHL.mat']],'CorrHL', 'methodID', 'meth','anovaTable_std_vs_N','anovaTable_N_vs_Nss','anovaTable_std_vs_Nnn','anovaTable_pred_vs_unpred');
saveas(s, [folderOut, [subj,'_SpatialCorrelation_BetasHL.svg']]);

