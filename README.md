# auditory_nordic

This repository accompanies a high-resolution 7T auditory dataset that is available on OPENNEURO: doi:10.18112/openneuro.ds004928.v1.0.0.
The data and code accompanies the manuscript now published in Imaging Neuroscience: https://doi.org/10.1162/imag_a_00270.

Note that the codes are written in Matlab and are rather poorly curated and specific for the data format we have used. We are working on better and more comprehensible versions of the code. 

Currently the input is data analyzed with BrainVoyager version 21.4 (Brain Inovation, Maastricht, The Netherlands). 

# Preprocessing scripts
In the folder NORDIC_preprocessing there are several scripts available. 
Dicoms are first converted to fmr (BV format) and then preprocessed using PreprocessingBV.m

In two participants,  we have done some extra distortion correction using ANTS because there were different distortions across runs (example scripts antsRegistration.sh and applyRegistration.sh) .
Moreover, to correct for geometric distortions using opposite phase encoding, we have used TOPUP from FSL (example scripts FMR2NII, topup.sh, applytopup.sh and NII2FMR).

Lastly, using a java script in BV, we create vtc files (BV format - example script VTCCreation.js).

# Processing after NORDIC - folder NORDIC_analyses

This folder contains several functions that are used in the main scripts.
For individual subjects, we have used explore_Nordic_effect_clean_v2.m.
For group analyses, we have used GroupAnalyses_NORDIC_v2.m
