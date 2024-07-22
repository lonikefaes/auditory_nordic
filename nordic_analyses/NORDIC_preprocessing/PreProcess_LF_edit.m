%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing script for fmr's in Brainvoyager 21.4 using Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

nr = 8;                                                             % Number of Runs
SubjNr = ['1'];

VerNOR = [{''}];                                                    % Enter Original, NORdef and NORnn as you please

nrVolumes = [230];

bv = actxserver('BrainVoyager21.BrainVoyagerScriptAccess.1');       % Specify which brainvoyager to use

for CurSub = 1:size(SubjNr,1)
    for j = 1:length(VerNOR)
        Subj = {['S',num2str(SubjNr)]};
        maindir = ['D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\S',num2str(SubjNr),'*_PC',VerNOR{j},'\Dicom\'];

        dirs = dir([maindir,'*RUN*_eja_ep2d_bold_pt8_TR1600_NoGap_NEW_012021']);        % Select proper directory
        
        for k = 1:length(dirs)   
            datadir2 = fullfile(dirs(k).folder,dirs(k).name);

            %% Cut fmrs to remove volumes acquired for NORDIC
            fmr_name = dir([datadir2, '\S', num2str(SubjNr),'*', VerNOR{j},'.fmr']);
            fmr = xff([datadir2, '\', fmr_name(1).name]);

            
            fmr.LoadSTC;
            newfmr = fmr.copyobject;
            sizefmr = size(fmr.Slice.STCData);
            newSize = sizefmr;
            newSize(3)= nrVolumes(1);                               % We only use volujme 1:'nrVolumes' to remove the noise volumes at the end
            newfmr.Slice.STCData = zeros(newSize);
            newfmr.NrOfVolumes = nrVolumes(1);
            newfmr.Slice.STCData = fmr.Slice.STCData(:,:,[1:nrVolumes(1)],:);
            
            newfmr.SaveAs([datadir2, '\S', num2str(SubjNr), '_MN_NG_PC_run',num2str(k), VerNOR{j}, '_Cut.fmr']);
            fmr.clearobject;
            clear fmr;
            newfmr.clearobject;
            clear newfmr;
            
            
            if strcmp(SubjNr(CurSub), '2')
                datadir_ref = fullfile(dirs(1).folder,dirs(4).name);
                fmrfile = dir(fullfile(datadir_ref,['*_PC_run4', VerNOR{j}, '.fmr']));              % Specify which run to do motion correction to - in our case mainly run 1, in two participants run 4
                mot_corr_ref = fullfile(fmrfile.folder,fmrfile.name);
            elseif strcmp(SubjNr(CurSub), '3')
                datadir_ref = fullfile(dirs(1).folder,dirs(4).name);
                fmrfile = dir(fullfile(datadir_ref,['*_PC_run4', VerNOR{j}, '.fmr'])); 
                mot_corr_ref = fullfile(fmrfile.folder,fmrfile.name);
%             elseif strcmp(SubjNr(CurSub), '1')
%                 datadir_ref = fullfile(dirs(k).folder,dirs(k).name);
%                 fmrfile = dir(fullfile(datadir_ref,['*_PC_run',num2str(k), VerNOR{j}, '.fmr']));
%                 mot_corr_ref = fullfile(fmrfile.folder,fmrfile.name);
            else
                datadir_ref = fullfile(dirs(1).folder,dirs(1).name);
                fmrfile = dir(fullfile(datadir_ref,['*_PC_run1', VerNOR{j}, '.fmr']));
                mot_corr_ref = fullfile(fmrfile.folder,fmrfile.name);
            end
            
            %% Start preprocessing
            %datadir = fullfile(dirs(CurSub).folder,dirs(CurSub).name);
            fmrfiles = dir(fullfile(datadir2,['\S', num2str(SubjNr),'_MN_NG_PC_',VerNOR{j}, '*Cut.fmr']));
            %         if i <=9
            [PathName FileName ext ] = fileparts(fullfile(fmrfiles.folder,fmrfiles.name)); % select the CG_OBJECTS file
            docFMR = bv.OpenDocument(fullfile(fmrfiles.folder,fmrfiles.name));
            
            
            % Set spatial and temporal parameters relevant for preprocessing
            if strcmp(SubjNr(CurSub), '5')                                              % Slightly different TR for NYU so for two subjects we use different TR
                if ~docFMR.TimeResolutionVerified
                    docFMR.TR = 1650;
                    docFMR.InterSliceTime = 39;
                    docFMR.TimeResolutionVerified = 1;
                end
                if ~docFMR.VoxelResolutionVerified
                    docFMR.PixelSizeOfSliceDimX = 0.801887;
                    docFMR.PixelSizeOfSliceDimY = 0.801887;
                    docFMR.SliceThickness = 0.8;
                    docFMR.GapThickness = 0;
                    docFMR.VoxelResolutionVerified = 1;
                end
            elseif strcmp(SubjNr(CurSub), '6')
                if ~docFMR.TimeResolutionVerified
                    docFMR.TR = 1650;
                    docFMR.InterSliceTime = 39;
                    docFMR.TimeResolutionVerified = 1;
                end
                if ~docFMR.VoxelResolutionVerified
                    docFMR.PixelSizeOfSliceDimX = 0.801887;
                    docFMR.PixelSizeOfSliceDimY = 0.801887;
                    docFMR.SliceThickness = 0.8;
                    docFMR.GapThickness = 0;
                    docFMR.VoxelResolutionVerified = 1;
                end
            else 
                if ~docFMR.TimeResolutionVerified
                    docFMR.TR = 1600;
                    docFMR.InterSliceTime = 38;
                    docFMR.TimeResolutionVerified = 1;
                end
                if ~docFMR.VoxelResolutionVerified
                    docFMR.PixelSizeOfSliceDimX = 0.804762;
                    docFMR.PixelSizeOfSliceDimY = 0.804762;
                    docFMR.SliceThickness = 0.8;
                    docFMR.GapThickness = 0;
                    docFMR.VoxelResolutionVerified = 1;
                end 
            end 
            % We save the new settings into the FMR file
            docFMR.Save;
            
            %%
            % Preprocessing step 2: Slice time correction
            if (docFMR.HasSliceTimeTable)
                docFMR.CorrectSliceTimingUsingTimeTable(2); % 2: Sinc Interpolation
            else
                docFMR.CorrectSliceTiming( 1, 0 );
                % First param: Scan order 0 -> Ascending, 1 -> Asc-Interleaved, 2 -> Asc-Int2,
                % 10 -> Descending, 11 -> Desc-Int, 12 -> Desc-Int2
                % Second param: Interpolation method: 0 -> trilinear, 1 -> sinc
            end
            ResultFileName = docFMR.FileNameOfPreprocessdFMR;
            docFMR.Close;
            docFMR = bv.OpenDocument( ResultFileName );
            % Preprocessing step 3: 3D motion correction
            %         if (nr==1)
            %             docFMR.CorrectMotion(2); %Trilinear Detection Sinc Interpolation (Sinc option in brain voyager)
            %             ResultFileName = docFMR.FileNameOfPreprocessdFMR;
            %         else
            docFMR.CorrectMotionTargetVolumeInOtherRunEx(mot_corr_ref,1,2,1,100,0,1);
            ResultFileName = docFMR.FileNameOfPreprocessdFMR;
            
            %         end
            %         1. target fmr file name
            %         2. target volume number usually the first volume, so 1.
            %         3. interpolation method use 1: trilinear
            %         4. use full data set if 1: yes, if 0: no
            %             5. maximum number of iterations this determines how often the estimates should be changed in order
            %             to find the values for the rotation and translation parameters; in the BrainVoyager QX user interface,
            %             this value is default set to 100.
            %             6. create movies these *.avi files cannot be generated yet, so the parameter can be set to 0.
            %             7. create extended log file this is very useful. if 1: yes, if 0: no
            % the current doc (input FMR) knows the name of the automatically saved output FMR
            %         docFMR.Remove; % close or remove input FMR
            % bvqx.PrintToLog(�Removed slice scan time corrected files instead of just closing...�)
            % docFMR.Close(); // close input FMR
            
            %% Spatial Smoothing - NOT USED 
            %         docFMR = bv.OpenDocument( ResultFileName );
            % % Open motion corrected file (output FMR) and assign to our doc var
            % % Preprocessing step 4: Spatial Gaussian Smoothing (not recommended
            % % for individual analysis with a 64x64 matrix)
            % docFMR.SpatialGaussianSmoothing( 4, �mm�); % FWHM value and unit
            % ResultFileName = docFMR.FileNameOfPreprocessdFMR;
            % docFMR.Close; % docFMR.Remove(); % close or remove input FMR
            % docFMR = bvqx.OpenDocument( ResultFileName );
            % Preprocessing step 5: Temporal High Pass Filter, includes Linear Trend
            % Removal
            
            %% Temporal Filtering
            docFMR.Close;
            docFMR = bv.OpenDocument( ResultFileName );
            docFMR.TemporalHighPassFilterGLMFourier(7); %2 sine cosine (DCT)
            ResultFileName = docFMR.FileNameOfPreprocessdFMR;
            docFMR.Close; % docFMR.Remove(); // close or remove input FMR
            
            %% Temporal Smoothing
            docFMR = bv.OpenDocument( ResultFileName );
            docFMR.TemporalGaussianSmoothing(2, 'TR');
            ResultFileName = docFMR.FileNameOfPreprocessdFMR;
            docFMR.Close; % docFMR.Remove(); // close or remove input FMR
            docFMR = bv.OpenDocument( ResultFileName );
            docFMR.Close;

        end
    end
end
