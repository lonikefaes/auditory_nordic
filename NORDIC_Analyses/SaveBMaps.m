function SaveBMaps (SaveBetaMapsMode,trial_x_run,connames,colnames,bBn,bAn,bAnNS,folderOut,roi_large,VTCInfo,subj)


if SaveBetaMapsMode == 1
    MapBn = [];
    MapAn = [];
    MapAn_NS = [];
    mapNames = {};
    count=1;
    for i = 1:max(trial_x_run)
        idx_t = find(trial_x_run==i);
        colnames_r = colnames(idx_t);
        for j=1:length(connames)
            select = contains(colnames_r, connames{j});
            idx_this = idx_t(select);
            MapBn = [MapBn,mean(bBn(idx_this,:))'];
            MapAn = [MapAn,mean(bAn(idx_this,:))'];
            MapAn_NS = [MapAn_NS,mean(bAnNS(idx_this,:))'];
            mapNames{count} = [connames{j},'_run',num2str(i)];
            count=count+1;
        end
    end
    create_vmp(folderOut,[subj, '_Betas_SingleRuns_bn.vmp'],VTCInfo,MapBn,roi_large,  mapNames);
    create_vmp(folderOut,[subj, '_Betas_SingleRuns_an.vmp'],VTCInfo,MapAn,roi_large,  mapNames);
    create_vmp(folderOut,[subj, '_Betas_SingleRuns_an_NoNoise.vmp'],VTCInfo,MapAn_NS,roi_large,  mapNames);
elseif SaveBetaMapsMode == 2
    for j=1:length(connames)
        select = contains(colnames, connames{j});
        create_vmp(folderOut,[subj, '_Betas_SingleTrial_',connames{j},'_bn.vmp'],VTCInfo,bBn(select,:)',roi_large,colnames(select));
        create_vmp(folderOut,[subj, '_Betas_SingleTrial_',connames{j},'_an.vmp'],VTCInfo,bAn(select,:)',roi_large,colnames(select));
        create_vmp(folderOut,[subj, '_Betas_SingleTrial_',connames{j},'_an_NoNoise.vmp'],VTCInfo,bAnNS(select,:)',roi_large,colnames(select));
    end
end
% Average Beta Maps
MapBn = [];
MapAn = [];
MapAn_NS = [];
mapNames = {};
count=1;

for j=1:length(connames)
    select = contains(colnames, connames{j});
    MapBn = [MapBn,mean(bBn(select,:))'];
    MapAn = [MapAn,mean(bAn(select,:))'];
    MapAn_NS = [MapAn_NS,mean(bAnNS(select,:))'];
    mapNames{count} = [connames{j}];
    count=count+1;
end

create_vmp(folderOut, [subj,'_BetasAverage_bn.vmp'],VTCInfo,MapBn,roi_large,  mapNames);
create_vmp(folderOut,[subj,'_Betas_Average_an.vmp'],VTCInfo,MapAn,roi_large,  mapNames);
create_vmp(folderOut,[subj, '_Betas_Average_an_NoNoise.vmp'],VTCInfo,MapAn_NS,roi_large,  mapNames);