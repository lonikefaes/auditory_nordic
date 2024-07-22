% Convert FMR to NII (Faruk's version)
clear all;

% path for the nifti files
pathOut = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\S11_MN_NG_PC_NoN\Dicom';
pathOut2 = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\S11_MN_NG_PC\Dicom\DistortionCorrection\';


% set fmr names
%fmrName{1}=[pathOut,'\DistortionCorrection\S11_MN_NG_PC_AP_3DMCTS.fmr'];
%fmrName{2}=[pathOut,'\DistortionCorrection\S11_MN_NG_PC_PA_3DMCTS.fmr'];
fmrName{3}=[pathOut,'\run1\S11_MN_NG_PC_run1_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr'];
fmrName{4}=[pathOut,'\run2\S11_MN_NG_PC_run2_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr'];
fmrName{5}=[pathOut,'\run3\S11_MN_NG_PC_run3_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr'];
fmrName{6}=[pathOut,'\run4\S11_MN_NG_PC_run4_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr'];
fmrName{7}=[pathOut,'\run5\S11_MN_NG_PC_run5_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr'];
fmrName{8}=[pathOut,'\run6\S11_MN_NG_PC_run6_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr'];
fmrName{9}=[pathOut,'\run7\S11_MN_NG_PC_run7_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr'];
fmrName{10}=[pathOut,'\run8\S11_MN_NG_PC_run8_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr'];

nrFmrs=length(fmrName);

% convert fmr to nii
disp('Converting fmr files to nii...')
for i=3:nrFmrs
    temp_fmr = xff(fmrName{i});
    
    % Output name parsing/creating
    [iPath, iName , iExt] = fileparts(fmrName{i});
    iExt = 'nii';
    outName = fullfile(pathOut2, [iName,'.',iExt]);
    
    temp_fmr.Write4DNifti(outName);
    temp_fmr.ClearObject;
    disp([outName,' ','created'])
end

disp('Done.')
