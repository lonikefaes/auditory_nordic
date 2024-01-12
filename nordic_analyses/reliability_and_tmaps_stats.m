function [] = reliability_and_tmaps_stats(rois,roiCount,roiNames,methCount,BetasCondrepSH, tCondrepSH , relBetas, relTs)

    tempData   = permute(BetasCondrepSH, [2,1,3] );
    tempMeth   = permute( methCount,     [2,1,3] ) ;    tempMeth  =    tempMeth(:,:);
    tempROI    = permute( roiCount,      [2,1,3] ) ;     tempROI  =    tempROI (:,:);

  
    t          = array2table(tempData(:,:));
    whithin    = table( categorical( tempMeth(1,:) )', categorical( tempROI(1,:))' ,'VariableNames', {'method', 'roi'} );
    rm         = fitrm(t,'Var1-Var15 ~ 1', 'WithinDesign', whithin); % Needs an interaction term, see below

    %% Contrasts:

     %         numeric 1*(         standard                   AND        ROI  =1) minus   convert to numeric by 1*(nordic AND ROI == 1)
    
    HG_Std_N   = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='1')+ 2.*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='1'));
    HG_N_Nnn   = find(1*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='1')+ 2.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='1')); 
    HG_Std_Nnn = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='1')+ 2.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='1')); 

    PP_Std_N   = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='2')- 1.*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='2'));
    PP_N_Nnn   = find(1*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='2')- 1.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='2')); 
    PP_Std_Nnn = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='2')- 1.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='2')); 

    PT_Std_N   = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='3')- 1.*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='3'));
    PT_N_Nnn   = find(1*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='3')- 1.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='3')); 
    PT_Std_Nnn = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='3')- 1.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='3')); 

    aSTG_Std_N   = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='4')- 1.*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='4'));
    aSTG_N_Nnn   = find(1*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='4')- 1.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='4')); 
    aSTG_Std_Nnn = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='4')- 1.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='4')); 

    pSTG_Std_N   = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='5')- 1.*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='5'));
    pSTG_N_Nnn   = find(1*(table2array(whithin(:,1))=='2' & table2array(whithin(:,2))=='5')- 1.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='5')); 
    pSTG_Std_Nnn = find(1*(table2array(whithin(:,1))=='1' & table2array(whithin(:,2))=='5')- 1.*(table2array(whithin(:,1))=='3' & table2array(whithin(:,2))=='5')); 


    %% Run statistics for betas

    % First we test for an interaction term
    tblBeta   = ranova(rm,'WithinModel',  'method*roi');        % If significant --> test all contrasts specified above

    BetaPerm = tempData(:,:)';
    Con = [HG_Std_N'; HG_N_Nnn'; HG_Std_Nnn'; PP_Std_N';PP_N_Nnn'; PP_Std_Nnn'; PT_Std_N'; PT_N_Nnn'; PT_Std_Nnn'; aSTG_Std_N'; aSTG_N_Nnn'; aSTG_Std_Nnn';...
        pSTG_Std_N'; pSTG_N_Nnn'; pSTG_Std_Nnn';];

    pBetas = Permutations_nonparametric(BetaPerm, size(BetaPerm,2), Con);

    pBetas_MC = pBetas*15;
%     aT_HG_Std_N   = ranova(rm,'WithinModel',  HG_Std_N);
%     aT_HG_N_Nnn  = ranova(rm,'WithinModel',  HG_N_Nnn);
%     aT_HG_Std_Nnn   = ranova(rm,'WithinModel', HG_Std_Nnn);
% 
%     aT_PP_Std_N   = ranova(rm,'WithinModel',  PP_Std_N);
%     aT_PP_N_Nnn  = ranova(rm,'WithinModel',  PP_N_Nnn);
%     aT_PP_Std_Nnn   = ranova(rm,'WithinModel', PP_Std_Nnn);
% 
%     aT_PT_Std_N   = ranova(rm,'WithinModel',  PT_Std_N);
%     aT_PT_N_Nnn  = ranova(rm,'WithinModel',  PT_N_Nnn);
%     aT_PT_Std_Nnn   = ranova(rm,'WithinModel', PT_Std_Nnn);
% 
%     aT_aSTG_Std_N   = ranova(rm,'WithinModel',  aSTG_Std_N);
%     aT_aSTG_N_Nnn  = ranova(rm,'WithinModel',  aSTG_N_Nnn);
%     aT_aSTG_Std_Nnn   = ranova(rm,'WithinModel', aSTG_Std_Nnn);
% 
%     aT_pSTG_Std_N   = ranova(rm,'WithinModel',  pSTG_Std_N);
%     aT_pSTG_N_Nnn  = ranova(rm,'WithinModel',  pSTG_N_Nnn);
%     aT_pSTG_Std_Nnn   = ranova(rm,'WithinModel', pSTG_Std_Nnn);
% 
%     pBetas = [aT_HG_Std_N{1,5} aT_HG_N_Nnn{1,5} aT_HG_Std_Nnn{1,5};...         % Put all p-values in a column and multiply by 15 for all the contrasts we test.
%         	  aT_PP_Std_N{1,5} aT_PP_N_Nnn{1,5} aT_PP_Std_Nnn{1,5};...
%               aT_PT_Std_N{1,5} aT_PT_N_Nnn{1,5} aT_PT_Std_Nnn{1,5};...
%               aT_aSTG_Std_N{1,5} aT_aSTG_N_Nnn{1,5} aT_aSTG_Std_Nnn{1,5};...
%               aT_pSTG_Std_N{1,5} aT_pSTG_N_Nnn{1,5} aT_pSTG_Std_Nnn{1,5}];
%     pBetas_MC = pBetas*15;

    %% Plot
    figure('Color',[1,1,1], 'Position',[281 19 1000 700])
    subplot(221)
    mb = boxchart(2.*roiCount(:),BetasCondrepSH(:),'GroupByColor', methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none');
    
    m    = 3;            % number of boxes per group
    xGap = 1 / m;

    m2 = xGap/2;

    ylim([-1 4]);
    h   = gca; set(h,'FontSize',12); 
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('$\beta$  \% change','interpreter','Latex')

    hold on

    plot(([2-xGap 2]), [3.5 3.5], '-k',  ([2-m2]), [3.6 3.6], '*k', 'MarkerSize', 3)
    %plot(([2-m2*0.5]), [3.6 3.6], '*k', 'MarkerSize', 3)
    plot(([2 2+xGap]), [3.7 3.7], '-k',  ([2+m2]), [3.8 3.8], '*k', 'MarkerSize', 3)
    %plot(([2+m2*0.5]), [3.8 3.8], '*k', 'MarkerSize', 3)
    plot(([2-xGap 2+xGap]), [3.9 3.9], '-k',  ([2]), [4 4], '*k', 'MarkerSize', 3)
    %plot(([2+m2*0.5]), [4 4], '*k', 'MarkerSize', 3)

    plot(([4-xGap 4]), [3.5 3.5], '-k',  ([4-m2]), [3.6 3.6], '*k', 'MarkerSize', 3)
    %plot(([4-m2*0.5]), [3.6 3.6], '*k', 'MarkerSize', 3)
    plot(([4 4+xGap]), [3.7 3.7], '-k',  ([4+m2]), [3.8 3.8], '*k', 'MarkerSize', 3)
    %plot(([4+m2*0.5]), [3.8 3.8], '*k', 'MarkerSize', 3)
    plot(([4-xGap 4+xGap]), [3.9 3.9], '-k',  ([4]), [4 4], '*k', 'MarkerSize', 3);
    %plot(([4+m2*0.5]), [4 4], '*k', 'MarkerSize', 3);


    plot(([6-xGap 6]), [3.5 3.5], '-k',  ([6-m2]), [3.6 3.6], '*k', 'MarkerSize', 3)
    %plot(([6-m2*0.5]), [3.6 3.6], '*k', 'MarkerSize', 3)
    plot(([6 6+xGap]), [3.7 3.7], '-k',  ([6+m2]), [3.8 3.8], '*k', 'MarkerSize', 3)
    %plot(([6+m2*0.5]), [3.8 3.8], '*k', 'MarkerSize', 3)
    plot(([6-xGap 6+xGap]), [3.9 3.9], '-k',  ([6]), [4 4], '*k', 'MarkerSize', 3)
    %plot(([6+m2*0.5]), [4 4], '*k', 'MarkerSize', 3)
    
    plot(([10-xGap 10]), [3.5 3.5], '-k',  ([10-m2]), [3.6 3.6], '*k', 'MarkerSize', 3)
    %plot(([10-m2*0.5]), [3.6 3.6], '*k', 'MarkerSize', 3)
    plot(([10 10+xGap]), [3.7 3.7], '-k',  ([10+m2]), [3.8 3.8], '*k', 'MarkerSize', 3)
    %plot(([10+m2*0.5]), [3.8 3.8], '*k', 'MarkerSize', 3)
    plot(([10-xGap 10+xGap]), [3.9 3.9], '-k',  ([10]), [4 4], '*k', 'MarkerSize', 3)
    %plot(([10+m2*0.5]), [4 4], '*k', 'MarkerSize', 3)


    %% T-values (plot B)
    tempData   = permute(tCondrepSH, [2,1,3] );    
    t          = array2table(tempData(:,:));
    rm          = fitrm(t,'Var1-Var15 ~ 1', 'WithinDesign' , whithin);

     % First we test for an interaction term
    tblT   = ranova(rm,'WithinModel',  'method*roi');        % If significant --> test all contrasts specified above

    TPerm = tempData(:,:)';
    
    [pT] = Permutations_nonparametric(TPerm, size(TPerm,2), Con);

    pT_MC = pT*15;

%     aT_HG_Std_N   = ranova(rm,'WithinModel',  HG_Std_N);
%     aT_HG_N_Nnn  = ranova(rm,'WithinModel',  HG_N_Nnn);
%     aT_HG_Std_Nnn   = ranova(rm,'WithinModel', HG_Std_Nnn);
% 
%     aT_PP_Std_N   = ranova(rm,'WithinModel',  PP_Std_N);
%     aT_PP_N_Nnn  = ranova(rm,'WithinModel',  PP_N_Nnn);
%     aT_PP_Std_Nnn   = ranova(rm,'WithinModel', PP_Std_Nnn);
% 
%     aT_PT_Std_N   = ranova(rm,'WithinModel',  PT_Std_N);
%     aT_PT_N_Nnn  = ranova(rm,'WithinModel',  PT_N_Nnn);
%     aT_PT_Std_Nnn   = ranova(rm,'WithinModel', PT_Std_Nnn);
% 
%     aT_aSTG_Std_N   = ranova(rm,'WithinModel',  aSTG_Std_N);
%     aT_aSTG_N_Nnn  = ranova(rm,'WithinModel',  aSTG_N_Nnn);
%     aT_aSTG_Std_Nnn   = ranova(rm,'WithinModel', aSTG_Std_Nnn);
% 
%     aT_pSTG_Std_N   = ranova(rm,'WithinModel',  pSTG_Std_N);
%     aT_pSTG_N_Nnn  = ranova(rm,'WithinModel',  pSTG_N_Nnn);
%     aT_pSTG_Std_Nnn   = ranova(rm,'WithinModel', pSTG_Std_Nnn);
%     
%     pT = [aT_HG_Std_N{1,5} aT_HG_N_Nnn{1,5} aT_HG_Std_Nnn{1,5};...         % Put all p-values in a column and multiply by 15 for all the contrasts we test.
%         	  aT_PP_Std_N{1,5} aT_PP_N_Nnn{1,5} aT_PP_Std_Nnn{1,5};...
%               aT_PT_Std_N{1,5} aT_PT_N_Nnn{1,5} aT_PT_Std_Nnn{1,5};...
%               aT_aSTG_Std_N{1,5} aT_aSTG_N_Nnn{1,5} aT_aSTG_Std_Nnn{1,5};...
%               aT_pSTG_Std_N{1,5} aT_pSTG_N_Nnn{1,5} aT_pSTG_Std_Nnn{1,5}];
% 
%     pT_MC = pT*15;

    %% Plot T's

    subplot(222);
    boxchart(2.*roiCount(:),tCondrepSH(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
    hold all;
    h = gca; set(h,'FontSize',12);
    ylim([prctile(tCondrepSH(:),1), prctile(tCondrepSH(:),99)])
    ylim([-1 4]);
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('t-statistic','interpreter','Latex')

    hold on

    plot(([2-xGap 2]), [3.5 3.5], '-k',  ([2-m2]), [3.6 3.6], '*k', 'MarkerSize', 3)
    %plot(([2-m2*0.5]), [3.6 3.6], '*k', 'MarkerSize', 3)
    plot(([2 2+xGap]), [3.7 3.7], '-k',  ([2+m2]), [3.8 3.8], '*k', 'MarkerSize', 3)
    plot(([2-xGap 2+xGap]), [3.9 3.9], '-k',  ([2]), [4 4], '*k', 'MarkerSize', 3)
    %plot(([2+m2*0.5]), [4 4], '*k', 'MarkerSize', 3)

    plot(([4-xGap 4+xGap]), [3.9 3.9], '-k',  ([4]), [4 4], '*k', 'MarkerSize', 3)

    plot(([6-xGap 6]), [3.5 3.5], '-k',  ([6-m2]), [3.6 3.6], '*k', 'MarkerSize', 3)
    plot(([6 6+xGap]), [3.7 3.7], '-k',  ([6+m2]), [3.8 3.8], '*k', 'MarkerSize', 3)
    plot(([6-xGap 6+xGap]), [3.9 3.9], '-k',  ([6]), [4 4], '*k', 'MarkerSize', 3)
    %plot(([6+m2*0.5]), [4 4], '*k', 'MarkerSize', 3)

    plot(([8-xGap 8]), [3.5 3.5], '-k',  ([8-m2]), [3.6 3.6], '*k', 'MarkerSize', 3)
    %plot(([8-m2*0.5]), [3.6 3.6], '*k', 'MarkerSize', 3)    
    plot(([8 8+xGap]), [3.7 3.7], '-k',  ([8+m2]), [3.8 3.8], '*k', 'MarkerSize', 3)
    plot(([8-xGap 8+xGap]), [3.9 3.9], '-k',  ([8]), [4 4], '*k', 'MarkerSize', 3)
   % plot(([8+m2*0.5]), [4 4], '*k', 'MarkerSize', 3)

    plot(([10-xGap 10]), [3.5 3.5], '-k',  ([10-m2]), [3.6 3.6], '*k', 'MarkerSize', 3)
    plot(([10 10+xGap]), [3.7 3.7], '-k',  ([10+m2]), [3.8 3.8], '*k', 'MarkerSize', 3)
    plot(([10-xGap 10+xGap]), [3.9 3.9], '-k',  ([10]), [4 4], '*k', 'MarkerSize', 3)
    %plot(([10+m2*0.5]), [4 4], '*k', 'MarkerSize', 3)
%     hleg = legend;
%     hleg.String(4:end)=[];
%     legend('Standard', 'NORDIC', 'NORDIC NN', 'Location', 'northeast', 'FontSize', 10, 'FontName', 'Arial')


    %% Reliability betas
   
   
    % bar( shBetas);title('$\beta$ split half corr (across vox)','interpreter','Latex'); hold all;
    % h = gca; set(h,'FontSize',12);set(h,'XTickLabel',connamesNoUnderscore,'TickLabelInterpreter','Latex')
    tempData   = permute(relBetas, [2,1,3] );    
    t          = array2table(tempData(:,:));
    rm          = fitrm(t,'Var1-Var15 ~ 1', 'WithinDesign' , whithin);
   
     % First we test for an interaction term
    tblRelB   = ranova(rm,'WithinModel',  'method*roi');        % If significant --> test all contrasts specified above

    RelBPerm = tempData(:,:)';
    
    [pRelB] = Permutations_nonparametric(RelBPerm, size(RelBPerm,2), Con);

     pRelB_MC = pRelB*15;

%     aT_HG_Std_N   = ranova(rm,'WithinModel',  HG_Std_N);
%     aT_HG_N_Nnn  = ranova(rm,'WithinModel',  HG_N_Nnn);
%     aT_HG_Std_Nnn   = ranova(rm,'WithinModel', HG_Std_Nnn);
% 
%     aT_PP_Std_N   = ranova(rm,'WithinModel',  PP_Std_N);
%     aT_PP_N_Nnn  = ranova(rm,'WithinModel',  PP_N_Nnn);
%     aT_PP_Std_Nnn   = ranova(rm,'WithinModel', PP_Std_Nnn);
% 
%     aT_PT_Std_N   = ranova(rm,'WithinModel',  PT_Std_N);
%     aT_PT_N_Nnn  = ranova(rm,'WithinModel',  PT_N_Nnn);
%     aT_PT_Std_Nnn   = ranova(rm,'WithinModel', PT_Std_Nnn);
% 
%     aT_aSTG_Std_N   = ranova(rm,'WithinModel',  aSTG_Std_N);
%     aT_aSTG_N_Nnn  = ranova(rm,'WithinModel',  aSTG_N_Nnn);
%     aT_aSTG_Std_Nnn   = ranova(rm,'WithinModel', aSTG_Std_Nnn);
% 
%     aT_pSTG_Std_N   = ranova(rm,'WithinModel',  pSTG_Std_N);
%     aT_pSTG_N_Nnn  = ranova(rm,'WithinModel',  pSTG_N_Nnn);
%     aT_pSTG_Std_Nnn   = ranova(rm,'WithinModel', pSTG_Std_Nnn);
%     
%     pRelB = [aT_HG_Std_N{1,5} aT_HG_N_Nnn{1,5} aT_HG_Std_Nnn{1,5};...         % Put all p-values in a column and multiply by 15 for all the contrasts we test.
%         	  aT_PP_Std_N{1,5} aT_PP_N_Nnn{1,5} aT_PP_Std_Nnn{1,5};...
%               aT_PT_Std_N{1,5} aT_PT_N_Nnn{1,5} aT_PT_Std_Nnn{1,5};...
%               aT_aSTG_Std_N{1,5} aT_aSTG_N_Nnn{1,5} aT_aSTG_Std_Nnn{1,5};...
%               aT_pSTG_Std_N{1,5} aT_pSTG_N_Nnn{1,5} aT_pSTG_Std_Nnn{1,5}];
% 
%     pRelB_MC = pRelB*15;

    subplot(223);
    boxchart(2.*roiCount(:),relBetas(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
    
    title('$\beta$ reliability','interpreter','Latex'); hold all;
    h = gca; set(h,'FontSize',12);
    ylim([-0.2 1])
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('split half corr','interpreter','Latex')
    hold on

  
    plot(([2-xGap 2]), [0.9 0.9], '-k',  ([2-m2]), [0.92 0.92], '*k', 'MarkerSize', 3)
    %plot(([2-m2*0.5]), [0.92 0.92], '*k', 'MarkerSize', 3)
    plot(([2-xGap 2+xGap]), [0.98 0.98], '-k',  ([2]), [1 1], '*k', 'MarkerSize', 3)
    %plot(([2+m2*0.5]), [1 1], '*k', 'MarkerSize', 3)

    plot(([6-xGap 6]), [0.9 0.9], '-k',  ([6-m2]), [0.92 0.92], '*k', 'MarkerSize', 3)
    %plot(([6-m2*0.5]), [0.92 0.92], '*k', 'MarkerSize', 3)
    plot(([6 6+xGap]), [0.94 0.94], '-k',  ([6+m2]), [0.96 0.96], '*k', 'MarkerSize', 3)
    plot(([6-xGap 6+xGap]), [0.98 0.98], '-k',  ([6]), [1 1], '*k', 'MarkerSize', 3)
    %plot(([6+m2*0.5]), [1 1], '*k', 'MarkerSize', 3)

    plot(([10-xGap 10+xGap]), [0.98 0.98], '-k',  ([10]), [1 1], '*k', 'MarkerSize', 3)


    %% Reliability of T
    tempData   = permute(relTs, [2,1,3] );    
    t          = array2table(tempData(:,:));
    rm          = fitrm(t,'Var1-Var15 ~ 1', 'WithinDesign' , whithin);

 % First we test for an interaction term
    tblRelT   = ranova(rm,'WithinModel',  'method*roi');        % If significant --> test all contrasts specified above

    RelTPerm = tempData(:,:)';
    
    [pRelT] = Permutations_nonparametric(RelTPerm, size(RelTPerm,2), Con);

    pRelT_MC = pRelT*15;

%     aT_HG_Std_N   = ranova(rm,'WithinModel',  HG_Std_N);
%     aT_HG_N_Nnn  = ranova(rm,'WithinModel',  HG_N_Nnn);
%     aT_HG_Std_Nnn   = ranova(rm,'WithinModel', HG_Std_Nnn);
% 
%     aT_PP_Std_N   = ranova(rm,'WithinModel',  PP_Std_N);
%     aT_PP_N_Nnn  = ranova(rm,'WithinModel',  PP_N_Nnn);
%     aT_PP_Std_Nnn   = ranova(rm,'WithinModel', PP_Std_Nnn);
% 
%     aT_PT_Std_N   = ranova(rm,'WithinModel',  PT_Std_N);
%     aT_PT_N_Nnn  = ranova(rm,'WithinModel',  PT_N_Nnn);
%     aT_PT_Std_Nnn   = ranova(rm,'WithinModel', PT_Std_Nnn);
% 
%     aT_aSTG_Std_N   = ranova(rm,'WithinModel',  aSTG_Std_N);
%     aT_aSTG_N_Nnn  = ranova(rm,'WithinModel',  aSTG_N_Nnn);
%     aT_aSTG_Std_Nnn   = ranova(rm,'WithinModel', aSTG_Std_Nnn);
% 
%     aT_pSTG_Std_N   = ranova(rm,'WithinModel',  pSTG_Std_N);
%     aT_pSTG_N_Nnn  = ranova(rm,'WithinModel',  pSTG_N_Nnn);
%     aT_pSTG_Std_Nnn   = ranova(rm,'WithinModel', pSTG_Std_Nnn);
%     
%     pRelT = [aT_HG_Std_N{1,5} aT_HG_N_Nnn{1,5} aT_HG_Std_Nnn{1,5};...         % Put all p-values in a column and multiply by 15 for all the contrasts we test.
%         	  aT_PP_Std_N{1,5} aT_PP_N_Nnn{1,5} aT_PP_Std_Nnn{1,5};...
%               aT_PT_Std_N{1,5} aT_PT_N_Nnn{1,5} aT_PT_Std_Nnn{1,5};...
%               aT_aSTG_Std_N{1,5} aT_aSTG_N_Nnn{1,5} aT_aSTG_Std_Nnn{1,5};...
%               aT_pSTG_Std_N{1,5} aT_pSTG_N_Nnn{1,5} aT_pSTG_Std_Nnn{1,5}];
%     pRelT_MC = pRelT*15;

    %% Plot reliability T
    subplot(224);
    boxchart(2.*roiCount(:),relTs(:),'GroupByColor',methCount(:),'WhiskerLineStyle' ,'none','MarkerStyle','none')
    
    title('t-value reliability','interpreter','Latex'); hold all;
    h = gca; set(h,'FontSize',12);
    ylim([-0.2 1])
    set(h,'XTick',unique(rois(:)))
    set(h,'XTickLabel',roiNames,'TickLabelInterpreter','Latex');
    ylabel('split half corr','interpreter','Latex')
    hold on

    plot(([2-xGap 2+xGap]), [0.98 0.98], '-k',  ([2]), [1 1], '*k', 'MarkerSize', 3)

    plot(([4-xGap 4]), [0.9 0.9], '-k',  ([4-m2]), [0.92 0.92], '*k', 'MarkerSize', 3)
    %plot(([4-m2*0.5]), [0.92 0.92], '*k', 'MarkerSize', 3)
    plot(([4 4+xGap]), [0.94 0.94], '-k',  ([4+m2]), [0.96 0.96], '*k', 'MarkerSize', 3)
    plot(([4-xGap 4+xGap]), [0.98 0.98], '-k',  ([4]), [1 1], '*k', 'MarkerSize', 3)
    %plot(([4+m2*0.5]), [1 1], '*k', 'MarkerSize', 3)

    plot(([6-xGap 6]), [0.9 0.9], '-k',  ([6-m2]), [0.92 0.92], '*k', 'MarkerSize', 3)
    plot(([6-xGap 6+xGap]), [0.98 0.98], '-k',  ([6]), [1 1], '*k', 'MarkerSize', 3)
    %plot(([6+m2*0.5]), [1 1], '*k', 'MarkerSize', 3)

    plot(([8-xGap 8]), [0.9 0.9], '-k',  ([8-m2]), [0.92 0.92], '*k', 'MarkerSize', 3)
    %plot(([8-m2*0.5]), [0.92 0.92], '*k', 'MarkerSize', 3)    
    plot(([8 8+xGap]), [0.94 0.94], '-k',  ([8+m2]), [0.96 0.96], '*k', 'MarkerSize', 3)
    %plot(([8+m2*0.5]), [0.96 0.96], '*k', 'MarkerSize', 3)
    plot(([8-xGap 8+xGap]), [0.98 0.98], '-k',  ([8]), [1 1], '*k', 'MarkerSize', 3)
    %plot(([8+m2*0.5]), [1 1], '*k', 'MarkerSize', 3)


    plot(([10-xGap 10]), [0.9 0.9], '-k',  ([10-m2]), [0.92 0.92], '*k', 'MarkerSize', 3)
    %plot(([10-m2*0.5]), [0.92 0.92], '*k', 'MarkerSize', 3)   
    plot(([10 10+xGap]), [0.94 0.94], '-k',  ([10+m2]), [0.96 0.96], '*k', 'MarkerSize', 3)
    %plot(([10+m2*0.5]), [0.96 0.96], '*k', 'MarkerSize', 3)
    plot(([10-xGap 10+xGap]), [0.98 0.98], '-k',  ([10]), [1 1], '*k', 'MarkerSize', 3)
    %plot(([10+m2*0.5]), [1 1], '*k', 'MarkerSize', 3)

   

