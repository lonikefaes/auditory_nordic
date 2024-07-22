#!/bin/bash

echo "This script will compute the topup"

echo "Merging files"
fslmerge -t up_down_phase S11_MN_NG_PC_AP_3DMCTS.nii S11_MN_NG_PC_PA_3DMCTS.nii
echo "Merging is done"

echo "Topup is running"
topup --imain=up_down_phase --datain=acqparams.txt --config=b0_faruk.cnf --out=topup_results
echo "Topup is computed"
