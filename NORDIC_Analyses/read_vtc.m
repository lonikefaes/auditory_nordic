function   [tsbn_t ,VTCInfo, tsnr ] = read_vtc(filename,roi)

    

    tsbn_tF   = xff( filename );
    VTCInfo   = extractVTCinfo( tsbn_tF);

    nrVoxRead = 1000;
    NrReadBlocks = floor(length(roi)/nrVoxRead) + 1;
    tsbn_t = zeros(VTCInfo.NrOfVolumes,length(roi));
    for j=1:NrReadBlocks
        
        if j<NrReadBlocks
            indexRead = roi((j-1)*nrVoxRead+1:j*nrVoxRead);
            tsbn_t(:,(j-1)*nrVoxRead+1:j*nrVoxRead) =  tsbn_tF.VTCData(:,indexRead);
        else
            if (j-1)*nrVoxRead+1<length(roi)
            indexRead = roi((j-1)*nrVoxRead+1:length(roi));
            else
            indexRead = [];
            end
            tsbn_t(:,(j-1)*nrVoxRead+1:length(roi)) =  tsbn_tF.VTCData(:,indexRead);
        end      
    end

%     tsbn_t    = tsbn_tF.VTCData;  
    tsnr      =  mean( tsbn_t )./std(tsbn_t);% compute tSNR per run
    % percent signal change
       tsbn_t    = 100.*tsbn_t./mean( tsbn_t );% this is standardization for percent signal change for each run
%    tsbn_t    = tsbn_t(:,roi);% this is standardization for percent signal change for each run
   tsbn_tF.clearobject;
   clear tsbn_tF;