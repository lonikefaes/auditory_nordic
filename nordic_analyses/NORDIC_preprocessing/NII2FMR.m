% Faruk's version
clc; clear all;

%%
% set path with nii files
pathIn  = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\S11_MN_NG_PC\Dicom\DistortionCorrection\';
pathOut = 'D:\LaminarfMRI_Audio\MN\NoGap\Post-Covid\S11_MN_NG_PC\Dicom\';

% set an fmr from the same subject to copy the header information
fmrName{1}=fullfile(pathOut, '\run1\S11_MN_NG_PC_run1_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr');
fmrName{2}=fullfile(pathOut, '\run2\S11_MN_NG_PC_run2_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr');
fmrName{3}=fullfile(pathOut, '\run3\S11_MN_NG_PC_run3_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr');
fmrName{4}=fullfile(pathOut, '\run4\S11_MN_NG_PC_run4_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr');
fmrName{5}=fullfile(pathOut, '\run5\S11_MN_NG_PC_run5_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr');
fmrName{6}=fullfile(pathOut, '\run6\S11_MN_NG_PC_run6_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr');
fmrName{7}=fullfile(pathOut, '\run7\S11_MN_NG_PC_run7_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr');
fmrName{8}=fullfile(pathOut, '\run8\S11_MN_NG_PC_run8_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.fmr');
nrFmrs=length(fmrName);

% set nii names created after spm or fsl manipulation
niiName{1}=fullfile(pathIn, 'S11_MN_NG_PC_run1_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii');
niiName{2}=fullfile(pathIn, 'S11_MN_NG_PC_run2_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii');
niiName{3}=fullfile(pathIn, 'S11_MN_NG_PC_run3_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii');
niiName{4}=fullfile(pathIn, 'S11_MN_NG_PC_run4_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii');
niiName{5}=fullfile(pathIn, 'S11_MN_NG_PC_run5_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii');
niiName{6}=fullfile(pathIn, 'S11_MN_NG_PC_run6_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii');
niiName{7}=fullfile(pathIn, 'S11_MN_NG_PC_run7_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii');
niiName{8}=fullfile(pathIn, 'S11_MN_NG_PC_run8_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii');

SeriesNr = [16 19 22 25 31 34 37 40];

% convert nii to fmr
nrFiles=length(fmrName);
disp('Converting nii files to fmr...')
for i=1:nrFiles
    % create temporary nii
    tempNii=xff(niiName{i});
    % Convert .nii back to _SCSTBL.fmr
    tempFmr=tempNii.Dyn3DToFMR;
    fmrProp=xff(fmrName{i});
    % Position information
    fmrProp.Slice.STCData=tempFmr.Slice.STCData;
    % Save
    [pathStrNii,fmrFileName,~] = fileparts(fmrName{i});
    fmrProp.SaveAs([pathOut,'\run', num2str(i), '\', fmrFileName(1:end),'_undist.fmr']);
    disp([pathStrNii, '\run', num2str(i), '\', fmrFileName(1:end), '_undist.fmr','created'])
end
disp('Done.')