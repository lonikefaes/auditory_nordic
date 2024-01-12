%%% Group analysis code for NORDIC versus Standard %%%

clc; clear all; close all;

%% Set main directories

datafolder = ['D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\'];

% outputfolder
folderOut = [datafolder,'ResultsGroupAnalysis',filesep];
if ~exist(folderOut)
    mkdir(folderOut);
end

allSubj = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11];

%% Load in all the data in the struct 'data'


for i = 1:length(allSubj)
    datadir = [datafolder, 'S', num2str(allSubj(i)), '_MN_NG_PC\ResultsScriptVolume\'];
    data{i}.CorrBetas       = load([datadir, 'S',num2str(allSubj(i)),'_MN_NG_PC_SpatialCorrelation_BetasAcrossCond.mat']);
    data{i}.RuntorunCorr    = load([datadir, 'S',num2str(allSubj(i)),'_MN_NG_PC_Run2runCorrelation_perMethod.mat']);
    data{i}.HL              = load([datadir, 'S',num2str(allSubj(i)),'_MN_NG_PC_SpatialCorrelation_BetasHL.mat']);
    %data{i}.RunCorr_AvBN    = load([datadir, 'S',num2str(allSubj(i)),'_MN_NG_PC_Run2runCorr_AverageBN.mat']);
    %data{i}.Tono            = load([datadir, 'S',num2str(allSubj(i)),'_MN_NG_PC_CorrelationTonotopy_run2Average.mat']);
end

%% Group average of spatial correlations --> Currently taking mean of correlations across runs and error is over participants
meanCorr = [];
for i = 1:length(allSubj)
    meanCorr(:,:,i) = mean(data{i}.CorrBetas.CorrBetas, [3]);
end

groupCorrBetas  = 0.*meanCorr; groupCorrBetas(1,:,:) = 1; groupCorrBetas(2,:,:) = 2;
groupCorrBetas2 = 0.*meanCorr; groupCorrBetas2(:,1,:) = 1; groupCorrBetas2(:,2,:) = 2; groupCorrBetas2(:,3,:) = 3; groupCorrBetas2(:,4,:) = 4;
groupCorrBetas2(:,5,:) = 5; groupCorrBetas2(:,6,:) = 6;
% check if the difference between conditions is statistically significant
%[aov1, tbl1, stats1] =  anova1(squeeze(  meanCorr(1,:,:) )' );
%[aov2, tbl2, stats2] =  anova1(squeeze(  meanCorr(2,:,:) )' );



%% Plot this figure
figure('Color',[1,1,1],'Position',[281 19 1160 500]),

subplot(1,14,1:6)

mb = boxchart(groupCorrBetas2(:),meanCorr(:), 'GroupByColor', groupCorrBetas(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');
ylabel('Group Average Spatial Correlation')
%xlabel('Condition','interpreter','Latex')
legend({'NORdef vs Original', 'NORnn vs Original'}, 'FontSize',10, 'FontName', 'Arial')
mb(1).BoxFaceColor = [144 103 167]./256;
mb(2).BoxFaceColor = [0.4660 0.6740 0.1880];
%title('Group Spatial Correlation')
xticks([1, 2, 3, 4, 5, 6])
xticklabels({'PredH', 'MispredH', 'UnpredH', 'PredL', 'MispredL', 'UnpredL'})
xlim([0 7])
hold on

% mkrs={'o','x','+', '*', '^', '.','v','<', '>', 'diamond'}; 
% %col = [0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840];
% col = {"#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#FF0000", "#0C5E1E", "#0C1F7C"};
% 
% set(groot, 'defaultTextHorizontalAlignment', 'center');
% m    = 2;            % number of boxes per group
% xGap = 1 / m;
% pp = [0.75 1.75 2.75 3.75 4.75 5.75];
% pt = [1.25 2.25 3.25 4.25 5.25 6.25];
% 
% 
% 
% for i = 1:length(meanCorr)
%     plot([pp], squeeze(meanCorr(1,:,i)), 'MarkerFaceColor', col{i}, 'MarkerEdgeColor', col{i}, 'Marker', mkrs(i), 'LineStyle','none')
%     plot([pt], squeeze(meanCorr(2,:,i)), 'MarkerFaceColor', col{i}, 'MarkerEdgeColor', col{i}, 'Marker', mkrs(i), 'LineStyle','none')
% end 

hleg = legend;
hleg.String(3:end)=[];
hleg.FontSize = 10;
hleg.FontName = 'Arial';

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';

% two way anova for comparing between conditions
%                         method                         x    conditions 

tempData = permute(meanCorr, [3, 2, 1]);
tempMeth   = permute(groupCorrBetas,     [3, 2, 1] ) ;   
tempMeth  =    tempMeth(:,:);
tempCond    = permute( groupCorrBetas2,      [3, 2, 1] ) ; 
tempCond = tempCond(:,:);
% 
% 
% t          = array2table(tempData(:,:)); % atanh fisher transform
% correlations to z-scale
t          = array2table(atanh(tempData(:,:)));
whithin    = table( categorical( tempMeth(1,:) )', categorical( tempCond(1,:))' ,'VariableNames', {'method', 'cond'} );
rm         = fitrm(t,'Var1-Var12 ~ 1', 'WithinDesign', whithin); % Needs an interaction term, see below

tbl   = ranova(rm,'WithinModel',  'method*cond');        % If significant --> test all contrasts specified above

%[H, pval2wayANOVA] = anovan(meanCorr(:), { categorical(groupCorrBetas(:)),     categorical(groupCorrBetas2(:))       });

% No signficant interaction, main effect of method, but not of cond.
% Average across condition and then test main effect of method.
av = squeeze(mean(tempData,2))';

[p] = Permutations_nonparametric(av, 10, [1 2]);



s = gcf;
% saveas(s,[folderOut, 'SpatialCorrelation_Group.fig']);
% save([folderOut,'SpatialCorrelation_Group.mat'],'meanCorr');
% saveas(s, [folderOut,'SpatialCorrelation_Group.svg']);

%% Average correlation run to run per method
methodID = zeros(3,length(allSubj));
for i = 1:length(allSubj)
    meanCorrR2R(:,i) = [mean(data{i}.RuntorunCorr.PPCor_BnH) mean(data{i}.RuntorunCorr.PPCor_AnH)  mean(data{i}.RuntorunCorr.PPCor_AnNSH)]';
    methodID(:,i)    = [1 2 3]';
end

per_method = methodID;

% figure('Color',[1,1,1],'Position',[281 19 1160 820]),

subplot(1,14,8:10)

% mb = boxchart( meanCorrR2R(:), 'GroupByColor', methodID(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none'); hold all;
mb = boxchart( meanCorrR2R(:),'GroupByColor', methodID(:),'WhiskerLineStyle' ,'none','MarkerStyle','none'); hold all;

ylabel('Correlation')
xlabel('Method')
%legend({'Standard', 'NORDIC', 'NORDIC NN'});
set(gca, 'XTick', [])
hold off
%title('Group Average Run to Run Correlation')

h = gcf;
hold on;

[A,B]      = sort(meanCorrR2R(1,:));

sorted_vec = [A; meanCorrR2R(2,B); meanCorrR2R(3,B)];
set(groot, 'defaultTextHorizontalAlignment', 'center');
m    = 3;            % number of boxes per group
xGap = 1 / m;
% plot([1-xGap 1 1+xGap], mean(meanCorrR2R,2)','o-')

for it = 1:length(allSubj)
    %     plot([mb(1).XData(it)  mb(2).XData(it) mb(3).XData(it)]', sorted_vec(:,it),'.','MarkerSize',20, 'Color', [0 0 0]+0.1*it); % no connecting lines
    %       plot( [1-xGap 1 1+xGap] ,sorted_vec(:,it),'.','MarkerSize',20, 'Color', [0 0 0]+0.1*it); % no connecting lines
    plot( [1-xGap 1 1+xGap] ,sorted_vec(:,it),'.','MarkerSize',10, 'Color', [1 1 1]-0.1*(0.9*it), 'MarkerSize',10 ); % no connecting lines

end

% hleg = legend;
% hleg.String(4:end)=[];
% xlim([1-2*xGap  1+2*xGap])
ylim([0 0.4])

hold on 
plot(([1-xGap 1+xGap]), [0.39 0.39], '-k',  ([1]), [0.4 0.4], '*k', 'MarkerSize', 6)


b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';


%[H, pval1wayANOVA, tbl] = anovan(meanCorrR2R(:), {  categorical(methodID(:))       });

[p] = Permutations_nonparametric(meanCorrR2R, length(meanCorrR2R),[1 2; 1 3]);
pval_CorrR2R = p * 2;


s = gcf;
% saveas(s,[folderOut, 'Run2runCorrelation_Betas.fig']);
% save([folderOut,'Run2runCorrelation_Betas.mat'],'meanCorrR2R');
% saveas(s, [folderOut,'Run2runCorrelation_Betas.svg']);
% anova for the difference
%[aov3, tbl3, stats3] = anova1(meanCorrR2R');


%% Average of LORO betas of CompH

methodID = zeros(3,length(allSubj));
meth = zeros(3,length(allSubj));
for i = 1:length(allSubj)
    %temp = [mean(data{i}.HL.CorrHL(:,1),[2])  mean(data{i}.HL.CorrHL(:,2,:),[3]) mean(data{i}.HL.CorrHL(:,3,:),[3])];
    %     temp = temp(:);
    temp = [data{i}.HL.CorrHL(1,:); data{i}.HL.CorrHL(3,:); data{i}.HL.CorrHL(5,:)];
    temp2 = mean(temp,[2]);
    meanHL(:,i) = temp2;
    methodID(:,i)    = [1 2 3]';    
    %methodIDa(:,i)   = [1 2 1 2 1 2]';% for anova use
    %meth(:,i)        = [1 2 3]';
end

%[aov4, tbl4, stats4] = anova1( vec( meanHL([1,3,5],:) ) , vec( meth([1,3,5],:) ) );     % comparing predictable High for each method
%[aov5, tbl5, stats5] = anova1( vec( meanHL([1,3,5]+1,:) ) , vec( meth([1,3,5]+1,:) ) ); % comparing predictable Low  for each method 


[A,B]      = sort(meanHL(1,:));
%[C,D]      = sort(meanHL(2,:));

sorted_vecH = [A; meanHL(2,B); meanHL(3,B)];
%sorted_vecL = [C; meanHL(4,D); meanHL(6,D)];

set(groot, 'defaultTextHorizontalAlignment', 'center');
m    = 3;            % number of boxes per group
xGap = 1 / m;

subplot(1,14,12:14)
mb = boxchart(meanHL(:), 'GroupByColor', methodID(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');
ylabel('Average Cross Validated Correlation')
xlabel('Method')
legend({'Original', 'NORdef', 'NORnn'});
set(gca, 'XTick', [])
ylim([0 0.6])

hold on

for it = 1:length(allSubj)
    %     plot([mb(1).XData(it)  mb(2).XData(it) mb(3).XData(it)]', sorted_vec(:,it),'.','MarkerSize',20, 'Color', [0 0 0]+0.1*it); % no connecting lines
    %       plot( [1-xGap 1 1+xGap] ,sorted_vec(:,it),'.','MarkerSize',20, 'Color', [0 0 0]+0.1*it); % no connecting lines
    %plot( [[1-xGap 2 2+xGap]], sorted_vecH(:,it),'.','MarkerSize',10, 'Color', [1 1 1]-0.1*(0.9*it), 'MarkerSize',10 ); % no connecting lines
    plot( [1-xGap 1 1+xGap], sorted_vecH(:,it),'.','MarkerSize',10, 'Color', [1 1 1]-0.1*(0.9*it), 'MarkerSize',10 ); % no connecting lines

    %plot( [1.66 4 6.33], sorted_vecL(:,it),'.','MarkerSize',10, 'Color', [1 1 1]-0.1*(0.9*it), 'MarkerSize',10 ); % no connecting lines


end



%xticks([1, 2, 3])
%xticklabels(['PredH'; 'PredH'; 'PredH'])

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';

hleg = legend;
hleg.String(4:end)=[];
hleg.FontSize = 10;
hleg.FontName = 'Arial';
hleg.Location = 'northeast';



%pval1wayANOVAa = anovan(meanHL(:), {  categorical(methodID(:)) });       %Before this was a two-way anova, but now we made it into one-way
[p] = Permutations_nonparametric(meanHL, length(meanHL),[1 2;1 3]);

pval_meanHL = p*2;

hold off
s = gcf;

% saveas(s,[folderOut, 'GroupFig_Correlations.fig']);
% saveas(s, [folderOut,'GroupFig_Correlations.svg']);
%saveas(s, [folderOut,'GroupFig_Correlations.png']);

%% Average of variance explained by design and what is left in residuals

Rbn =[];
Ran = [];
RanNS = [];
Rbres_refSSY = [];
RbresNS_refSSY = [];
methodIDres = [];
for i = 1:length(allSubj)
    dataHist = [datafolder, 'S', num2str(allSubj(i)), '_MN_NG_PC\ResultsScriptVolume\'];
    data     = load([dataHist, 'S',num2str(allSubj(i)),'_MN_NG_PC_SumOfS_Model.mat']);

    Rbn   =   [Rbn;   (data.Rbn')];
    Ran   =   [Ran;  (data.Ran')];
    RanNS =   [RanNS;(data.RanNS')];
    Rbres_refSSY =   [Rbres_refSSY; (data.Rbres_refSSY')];
    RbresNS_refSSY = [RbresNS_refSSY; (data.RbresNS_refSSY')];

    EVdesign(:,i) = [mean(data.Rbn')  mean(data.Ran') mean(data.RanNS')]';
    EVres(:,i)    = [mean(data.Rbres_refSSY') mean(data.RbresNS_refSSY') ]';
    methodIDres(:,i) = [1 2]';

end

%[aov4, tbl4, stats4] =  anova1(EVdesign');  %ANOVA for variance explained by design
%multcompare(stats4)                         %Pairwise comparisons
%[aov5, tbl5, stats5] =  anova1(EVdesign');
[p] = Permutations_nonparametric(EVdesign, length(EVdesign),[1 2; 1 3; 2 3]);     % Do we do permutations?

p_EVdesign = p*3;

%pval1wayANOVA = anovan(EVdesign(:), {  categorical(methodID(:))       });


%x = EVres(1,:);     % NORDIC residuals for all 10 participants.
%y = EVres(2,:);     % NORDIC No NOISE residuals for all 10 participants
%[H, P, CI, stats] = ttest(x,y)              % t-test for checking what is removed from the data??
[p_EVres] = Permutations_nonparametric(EVres, length(EVdesign),[1 2]);     % Do we do permutations?


figure('Color',[1,1,1],'Position',[281 19 1100 600]),
subplot(121);
mb = boxchart(EVdesign(:), 'GroupByColor', per_method(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');
ylabel('$\frac{Design\ SS}{Total\ SS}$', 'interpreter','Latex','FontSize',20)
xlabel('Method');
legend({'Original', 'NORdef', 'NORnn'});
set(gca, 'XTick', [])

h = gcf;
hold on;

[A,B]      = sort(EVdesign(1,:));

sorted_vec = [A; EVdesign(2,B); EVdesign(3,B)];
set(groot, 'defaultTextHorizontalAlignment', 'center');
m    = 3;            % number of boxes per group
xGap = 1 / m;

m2 = xGap/3;


for it = 1:length(allSubj)
    plot( [1-xGap 1 1+xGap] ,sorted_vec(:,it),'.','MarkerSize',20, 'Color', [1 1 1]-0.1*(0.9*it), 'MarkerSize',30 ); % no connecting lines
end

ylim([0.3 0.48])

plot(([1-xGap 1]), [0.445 0.445], '-k',  ([1-m2]), [0.45 0.45], '*k', 'MarkerSize', 6)
plot(([1-xGap 1]), [0.445 0.445], '-k',  ([1-2*m2]), [0.45 0.45], '*k', 'MarkerSize', 6)

plot(([1 1+xGap]), [0.455 0.455], '-k',  ([1+m2]), [0.46 0.46], '*k', 'MarkerSize', 6)
plot(([1 1+xGap]), [0.455 0.455], '-k',  ([1+2*m2]), [0.46 0.46], '*k', 'MarkerSize', 6)

plot(([1-xGap 1+xGap]), [0.465 0.465], '-k',  ([1-m2/2]), [0.47 0.47], '*k', 'MarkerSize', 6)
plot(([1-xGap 1+xGap]), [0.465 0.465], '-k',  ([1+m2/2]), [0.47 0.47], '*k', 'MarkerSize', 6)


hleg = legend;
hleg.String(4:end)=[];
hleg.FontSize = 10;
hleg.FontName = 'Arial';

% b = gca;
%b.XAxis.FontSize = 12;
%b.XAxis.FontName = 'Arial';
%b.YAxis.FontSize = 12;
%b.YAxis.FontName = 'Arial';

subplot(122);
mb = boxchart(EVres(:), 'GroupByColor', methodIDres(:), 'WhiskerLineStyle' ,'none','MarkerStyle','none');
mb(1).BoxFaceColor = [144 103 167]./256;
mb(2).BoxFaceColor = [0.4660 0.6740 0.1880];
ylabel('$\frac{Design\ SS\ Residuals}{Total\ SS\ Original}$', 'interpreter','Latex','fontsize', 20)
xlabel('Method')
legend({'NORdef vs Original', 'NORnn vs Original'});
set(gca, 'XTick', [])

h = gcf;
hold on;

[A,B]      = sort(EVres(1,:));

sorted_vec = [A; EVres(2,B)];
set(groot, 'defaultTextHorizontalAlignment', 'center');
m    = 4;            % number of boxes per group
xGap = 1 / m;
m2 = xGap/4;


for it = 1:length(allSubj)
    plot( [1-xGap 1+xGap] ,sorted_vec(:,it),'.','MarkerSize',20, 'Color', [1 1 1]-0.1*(0.9*it), 'MarkerSize',30 ); % no connecting lines

end

ylim([0.06 0.24])

hold on 

plot(([1-xGap 1+xGap]), [0.23 0.23], '-k',  ([1-m2]), [0.24 0.24], '*k', 'MarkerSize', 6)
plot(([1+m2]), [0.24 0.24], '*k', 'MarkerSize', 6)


hleg = legend;
hleg.String(3:end)=[];
hleg.FontSize = 10;
hleg.FontName = 'Arial';

b = gca;
%b.XAxis.FontSize = 12;
%b.XAxis.FontName = 'Arial';
% b.YAxis.FontSize = 12;
% b.YAxis.FontName = 'Arial';

h = gcf;
% saveas(h,[folderOut, 'Average_SumofS_Model.fig']);
% saveas(h, [folderOut, 'Average_SumofS_Model.svg']);


%% Group average of betas, t-values and consistency
for i = 1:length(allSubj)
    datadir = [datafolder, 'S', num2str(allSubj(i)), '_MN_NG_PC\ResultsScriptVolume\'];
    data     = load([datadir, 'S',num2str(allSubj(i)),'_MN_NG_PC_BetasANDt_changeANDReliability.mat']);
    BetasCondrepSH(:,i,:) = squeeze(mean(data.BetasCondrepSH,2));
    tCondrepSH(:,i,:)     = squeeze(mean(data.tCondrepSH,2));
    relBetas(:,i,:)       = squeeze(mean(data.relBetas,2));
    relTs(:,i,:)          = squeeze(mean(data.relTs,2));
    roiCount(:,i,:)       = repmat([1:5]',1,3);
    methCount(:,i,:)      = repmat([1,2,3],5,1 );
end
rois        = 2*(1:5);
roiNames    = {'HG','PP','PT','aSTG','pSTG',};
reliability_and_tmaps_stats(rois,roiCount,roiNames,methCount,BetasCondrepSH, tCondrepSH , relBetas, relTs);

h = gcf;

b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';

% saveas(h,[folderOut, 'Average_BetaAndT_Reliability.fig']);
% saveas(h, [folderOut, 'Average_BetaAndT_Reliability.svg']);
% saveas(h, [folderOut, 'Average_BetaAndT_Reliability.png']);
