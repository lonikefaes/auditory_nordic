% General ANTs
##!/bin/bash

# This is an example script in which we apply ANTS to do non-linear distortion correction.
# We used this because two participants had distortions across runs (in 2 or 3 runs) that we not related the usual distortions for which we use TOPUP.
#
  mydata=/mnt/d/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S7_MN_NG_PC/run8/


  myTarget=${mydata}/run1_cut.nii
  mySource=${mydata}/run8_cut.nii
  myITK=${mydata}
  outfld=${mydata}
  cd ${outfld}
  cp ${mySource} ./Moving_run8.nii
  cp ${myTarget} ./Target_run8.nii

  #cp ${myITK}/initial_matrix_ITK.txt .


  # Coregistration done in 2 steps
  echo "*****************************************"
  echo "************* starting with ANTS ********"
  echo "*****************************************"
  antsRegistration \
   --verbose 1 \
   --dimensionality 3 \
   --float 1 \
   --output [run8_registered_,run8_registered_Warped.nii.gz,run8_registered_InverseWarped.nii.gz] \
   --interpolation BSpline[5] \
   --use-histogram-matching 0 \
   --winsorize-image-intensities [0.005,0.995] \
   --transform Rigid[0.05] \
   --metric MI[Target_run8.nii,Moving_run8.nii,0.7,32,Regular,0.1] \
   --convergence [1000x500,1e-6,10] \
   --shrink-factors 4x2 \
   --smoothing-sigmas 1x0vox \
   --transform Affine[0.1] \
   --metric MI[Target_run8.nii,Moving_run8.nii,0.7,32,Regular,0.1] \
   --convergence [1000x500,1e-6,10] \
   --shrink-factors 4x2 \
   --smoothing-sigmas 1x0vox \
   --transform SyN[0.1,2,0] \
   --metric CC[Target_run8.nii,Moving_run8.nii,1,2] \
   --convergence [500x100,1e-6,10] \
   --shrink-factors 2x1 \
   --smoothing-sigmas 1x0vox
