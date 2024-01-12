function SaveTonoRes(thF,connames,colnames,bBn,bAn,Fbn,folderOut,roi_large,VTCInfo,subj)

% no selection - all voxels in large ROI
filename = ['_Tono_noThRes.vmp'];
Map = [];
count = 1;
mapNames = [];
for i=1:size(connames,1)
    tonoBn = [];
    tonoAn = [];
    temp = connames{i,1};
    condName = temp(1:end-2);
    if i==1
        mapNames = {['Tono_',condName,'_res'],['Tono_',condName,'_resNN']};
    else
        mapNames = [mapNames,['Tono_',condName,'_res'],['Tono_',condName,'_resNN']];
    end
    for j=1:size(connames,2)
        select = contains(colnames, connames{i,j});
        tonoBn(j,:) = mean(bBn(select,:));
        tonoAn(j,:) = mean(bAn(select,:));
    end
    tonoBn = zscore(tonoBn,[],2);
    tonoAn = zscore(tonoAn,[],2);
   
    [~,indM_bn] = max(tonoBn);
    [~,indM_an] = max(tonoAn);
    
    indM_bn(indM_bn==1) = -10; % high
    indM_bn(indM_bn==2) = 10;  % low
    indM_an(indM_an==1) = -10; % high
    indM_an(indM_an==2) = 10;  % low
   
    Map(:,count) =indM_bn;
    Map(:,count+1) =indM_an;
   
    count = count+2;
end
 create_vmp_tono(folderOut,[subj, filename],VTCInfo,Map,roi_large, mapNames);

% select on bn - F>thF
indexVox = find(Fbn>thF);
filename = ['_Tono_ThOnRes.vmp'];
Map = [];
count = 1;
mapNames = [];
for i=1:size(connames,1)
    tonoBn = [];
    tonoAn = [];
   
    temp = connames{i,1};
    condName = temp(1:end-2);
    if i==1
        mapNames = {['Tono_',condName,'_res'],['Tono_',condName,'_resNN']};
    else
        mapNames = [mapNames,['Tono_',condName,'_res'],['Tono_',condName,'_resNN']];
    end
    for j=1:size(connames,2)
        select = contains(colnames, connames{i,j});
        tonoBn(j,:) = mean(bBn(select,indexVox));
        tonoAn(j,:) = mean(bAn(select,indexVox));
    end
    tonoBn = zscore(tonoBn,[],2);
    tonoAn = zscore(tonoAn,[],2);
   
    [~,indM_bn] = max(tonoBn);
    [~,indM_an] = max(tonoAn);
 
    indM_bn(indM_bn==1) = -10; % high
    indM_bn(indM_bn==2) = 10;  % low
    indM_an(indM_an==1) = -10; % high
    indM_an(indM_an==2) = 10;  % low
   
    Map(:,count) =indM_bn;
    Map(:,count+1) =indM_an;
   
    count = count+2;
end
 create_vmp_tono(folderOut,[subj,filename],VTCInfo,Map,roi_large(indexVox), mapNames);


