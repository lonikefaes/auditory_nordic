% General ANTs
#!/bin/bash

set -m

echo "I am now applying the transformation to run 8 original"

mydata=/mnt/d/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S7_MN_NG_PC/run8/run8/run8


myTarget=${mydata}/run1_cut.nii
#mySource=${mydata}/run6_cut.nii
myITK=${mydata}
outfld=${mydata}
cd ${outfld}
antsApplyTransforms -d 3 -e 3 -i S7_MN_NG_PC_run8_Cut_SCSTBL_3DMAS_3DMCTS.nii -o S7_MN_NG_PC_run8_Cut_SCSTBL_3DMAS_3DMCTS_warped.nii -r run1_cut.nii -t run8_registered_1Warp.nii.gz -t run8_registered_0GenericAffine.mat

