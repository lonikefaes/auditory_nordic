function SaveTono(thF,connames,colnames,bBn,bAn,bAnNS,Fbn,Fan,FanNS,folderOut,roi_large,VTCInfo, subj)

% no selection - all voxels in large ROI
filename = ['_Tono_noTh.vmp'];
Map = [];
count = 1;
mapNames = [];
for i=1:size(connames,1)
    tonoBn = [];
    tonoAn = [];
    tonoAnNS = [];
    temp = connames{i,1};
    condName = temp(1:end-2);
    if i==1
        mapNames = {['Tono_',condName,'_bn'],['Tono_',condName,'_an'],['Tono_',condName,'_an_NoNoise']};
    else
        mapNames = [mapNames,['Tono_',condName,'_bn'],['Tono_',condName,'_an'],['Tono_',condName,'_an_NoNoise']];
    end
    for j=1:size(connames,2)
        select = contains(colnames, connames{i,j});
        tonoBn(j,:) = mean(bBn(select,:));
        tonoAn(j,:) = mean(bAn(select,:));
        tonoAnNS(j,:) = mean(bAnNS(select,:));
    end
    

    tonoBn = zscore(tonoBn,[],2);
    tonoAn = zscore(tonoAn,[],2);
    tonoAnNS = zscore(tonoAnNS,[],2);
    [~,indM_bn] = max(tonoBn);
    [~,indM_an] = max(tonoAn);
    [~,indM_anNS] = max(tonoAnNS);
    indM_bn(indM_bn==1) = -10; % high
    indM_bn(indM_bn==2) = 10;  % low
    indM_an(indM_an==1) = -10; % high
    indM_an(indM_an==2) = 10;  % low
    indM_anNS(indM_anNS==1) = -10; % high
    indM_anNS(indM_anNS==2) = 10;  % low
    Map(:,count) =indM_bn;
    Map(:,count+1) =indM_an;
    Map(:,count+2) =indM_anNS;
    count = count+3;
end
 create_vmp_tono(folderOut,[subj, filename],VTCInfo,Map,roi_large, mapNames);

% select on bn - F>thF
indexVox = find(Fbn>thF);
filename = ['_Tono_ThOnBn.vmp'];
Map = [];
count = 1;
mapNames = [];
for i=1:size(connames,1)
    tonoBn = [];
    tonoAn = [];
    tonoAnNS = [];
    temp = connames{i,1};
    condName = temp(1:end-2);
    if i==1
        mapNames = {['Tono_',condName,'_bn'],['Tono_',condName,'_an'],['Tono_',condName,'_an_NoNoise']};
    else
        mapNames = [mapNames,['Tono_',condName,'_bn'],['Tono_',condName,'_an'],['Tono_',condName,'_an_NoNoise']];
    end
    for j=1:size(connames,2)
        select = contains(colnames, connames{i,j});
        tonoBn(j,:) = mean(bBn(select,indexVox));
        tonoAn(j,:) = mean(bAn(select,indexVox));
        tonoAnNS(j,:) = mean(bAnNS(select,indexVox));
    end
    tonoBn = zscore(tonoBn,[],2);
    tonoAn = zscore(tonoAn,[],2);
    tonoAnNS = zscore(tonoAnNS,[],2);
    [~,indM_bn] = max(tonoBn);
    [~,indM_an] = max(tonoAn);
    [~,indM_anNS] = max(tonoAnNS);
    indM_bn(indM_bn==1) = -10; % high
    indM_bn(indM_bn==2) = 10;  % low
    indM_an(indM_an==1) = -10; % high
    indM_an(indM_an==2) = 10;  % low
    indM_anNS(indM_anNS==1) = -10; % high
    indM_anNS(indM_anNS==2) = 10;  % low
    Map(:,count) =indM_bn;
    Map(:,count+1) =indM_an;
    Map(:,count+2) =indM_anNS;
    count = count+3;
end
 create_vmp_tono(folderOut,[subj, filename],VTCInfo,Map,roi_large(indexVox), mapNames);




% select on F>thF - each processing type independently
indexVox1 = find(Fbn>thF);
indexVox2 = find(Fan>thF);
indexVox3 = find(FanNS>thF);
filename = ['_Tono_ThInd.vmp'];
Map = zeros(length(roi_large),9);
count = 1;
mapNames = [];
for i=1:size(connames,1)
    tonoBn = [];
    tonoAn = [];
    tonoAnNS = [];
    temp = connames{i,1};
    condName = temp(1:end-2);
    if i==1
        mapNames = {['Tono_',condName,'_bn'],['Tono_',condName,'_an'],['Tono_',condName,'_an_NoNoise']};
    else
        mapNames = [mapNames,['Tono_',condName,'_bn'],['Tono_',condName,'_an'],['Tono_',condName,'_an_NoNoise']];
    end
    for j=1:size(connames,2)
        select = contains(colnames, connames{i,j});
        tonoBn(j,:) = mean(bBn(select,indexVox1));
        tonoAn(j,:) = mean(bAn(select,indexVox2));
        tonoAnNS(j,:) = mean(bAnNS(select,indexVox3));
    end
    tonoBn = zscore(tonoBn,[],2);
    tonoAn = zscore(tonoAn,[],2);
    tonoAnNS = zscore(tonoAnNS,[],2);
    [~,indM_bn] = max(tonoBn);
    [~,indM_an] = max(tonoAn);
    [~,indM_anNS] = max(tonoAnNS);
    indM_bn(indM_bn==1) = -10; % high
    indM_bn(indM_bn==2) = 10;  % low
    indM_an(indM_an==1) = -10; % high
    indM_an(indM_an==2) = 10;  % low
    indM_anNS(indM_anNS==1) = -10; % high
    indM_anNS(indM_anNS==2) = 10;  % low
    Map(indexVox1,count) =indM_bn;
    Map(indexVox2,count+1) =indM_an;
    Map(indexVox3,count+2) =indM_anNS;
    count = count+3;
end
 create_vmp_tono(folderOut,[subj,filename],VTCInfo,Map,roi_large, mapNames);