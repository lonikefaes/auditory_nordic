#!/bin/bash

cd /mnt/d/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/DistortionCorrection/

run_num=('1 2 3 4 5 6 7 8')
for j in ${run_num[@]}; do

  echo "Changing NaNs to zeros"

  fslmaths 'S11_MN_NG_PC_run'$j'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.nii' -nan 'S11_MN_NG_PC_run'$j'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.nii.gz'
  rm 'S11_MN_NG_PC_run'$j'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.nii'

  echo "Applying topup"

  applytopup -i 'S11_MN_NG_PC_run'$j'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.nii.gz' -a acqparams.txt --topup=topup_results --inindex=1 --method=jac --interp=spline --verbose --out='S11_MN_NG_PC_run'$j'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii.gz'

  gunzip 'S11_MN_NG_PC_run'$j'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii.gz'
  rm -rf gunzip 'S11_MN_NG_PC_run'$j'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.nii.gz'

  rm 'S11_MN_NG_PC_run'$j'_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp.nii.gz'


done
