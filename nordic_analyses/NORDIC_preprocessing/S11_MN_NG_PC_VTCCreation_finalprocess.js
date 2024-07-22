// Batch script for VTC creation in BV 21.4

// From native functional data space to an anatomical data to ACPC

// For more information, check the scripting guide:

// http://support.brainvoyager.com/automation-aamp-development/46-writing-scripts/133-scripting-reference-guide.html



var bvqx =BrainVoyager;

bvqx.PrintToLog("---");



//----------------------------------------------------------------

// Enter the following values:

var OutDataPath = "D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/";                     // Path to where output files are saved

var NumberOfRuns = 24;                                                                          // 8 runs for this participant, and 3 methods (Original, NORdef, NORnn)

var AnatomicalData = "D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/MR-SE004-t1_mp2rage_sag_p3_0.75mm_UNI_Images/S11_MN_NG_PC_UNI_ISO-0.7_IIHC_pt4_ACPC.vmr";     // Anatomical file in ACPC space



var IA_TRFList =  OutDataPath +"Dicom/run1/S11_MN_NG_PC_run1_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist-TO-S11_MN_NG_PC_UNI_ISO-0.7_IIHC_pt4_IA.trf";                     // Initial alignment file.

var FA_TRFList =  OutDataPath + "Dicom/run1/S11_MN_NG_PC_run1_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist-TO-S11_MN_NG_PC_UNI_ISO-0.7_IIHC_pt4_BBR_FA.trf";                // Fine alignment file.

var ACPC_TRFList = OutDataPath + "Dicom/MR-SE004-t1_mp2rage_sag_p3_0.75mm_UNI_Images/S11_MN_NG_PC_UNI_ISO-0.7_IIHC_pt4_ACPC.trf";                                               // Transformation file for anatomy from NATIVE to ACPC space.





// fmr file names usually get quite long after preprocessing

var InFmrList = [

"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/run1/S11_MN_NG_PC_run1_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/run2/S11_MN_NG_PC_run2_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/run3/S11_MN_NG_PC_run3_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/run4/S11_MN_NG_PC_run4_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/run5/S11_MN_NG_PC_run5_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/run6/S11_MN_NG_PC_run6_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/run7/S11_MN_NG_PC_run7_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC/Dicom/run8/S11_MN_NG_PC_run8_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NORDIC/Dicom/run1/S11_MN_NG_PC_run1_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NORDIC/Dicom/run2/S11_MN_NG_PC_run2_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NORDIC/Dicom/run3/S11_MN_NG_PC_run3_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NORDIC/Dicom/run4/S11_MN_NG_PC_run4_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NORDIC/Dicom/run5/S11_MN_NG_PC_run5_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NORDIC/Dicom/run6/S11_MN_NG_PC_run6_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NORDIC/Dicom/run7/S11_MN_NG_PC_run7_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NORDIC/Dicom/run8/S11_MN_NG_PC_run8_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NoN/Dicom/run1/S11_MN_NG_PC_run1_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NoN/Dicom/run2/S11_MN_NG_PC_run2_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NoN/Dicom/run3/S11_MN_NG_PC_run3_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NoN/Dicom/run4/S11_MN_NG_PC_run4_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NoN/Dicom/run5/S11_MN_NG_PC_run5_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NoN/Dicom/run6/S11_MN_NG_PC_run6_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NoN/Dicom/run7/S11_MN_NG_PC_run7_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr",
"D:/LaminarfMRI_Audio/MN/NoGap/Post-Covid/S11_MN_NG_PC_NoN/Dicom/run8/S11_MN_NG_PC_run8_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.fmr"];                             // All fmr files that need to be converted to VTC. 
                                                                                                                                                                                // Preprocessing included cutting last volumes, slice scan time correction, moco, temporal filtering and smoothing of 2 datapoints, and distortion correction.






// Output Talairach vtc file names

var OutVTCNameList = [

"S11_MN_NG_PC_run1_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run2_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run3_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run4_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run5_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run6_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run7_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run8_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run1_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run2_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run3_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run4_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run5_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run6_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run7_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run8_NORDIC_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run1_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run2_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run3_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run4_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run5_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run6_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run7_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc",
"S11_MN_NG_PC_run8_NoN_Cut_SCSTBL_3DMCTS_LTR_THPGLMF7c_TDTS2.0dp_undist.vtc"];                              // VTC filenames


//---------------------------------------------------------


for(count = 1; count <= NumberOfRuns; count++){


// Open an anatomical data

var docVMR = bvqx.OpenDocument(AnatomicalData);



// Specify bounding box for target VTC

docVMR.ExtendedTALSpaceForVTCCreation = false;  // optional

docVMR.UseBoundingBoxForVTCCreation = true; // optional

docVMR.TargetVTCBoundingBoxXStart = 170;

docVMR.TargetVTCBoundingBoxXEnd   = 590;

docVMR.TargetVTCBoundingBoxYStart = 180;

docVMR.TargetVTCBoundingBoxYEnd   = 640;

docVMR.TargetVTCBoundingBoxZStart = 250;

docVMR.TargetVTCBoundingBoxZEnd   = 500;

    var InFmrPath = InFmrList[count-1];

    var OutVTCName = OutDataPath + "VTC/" + OutVTCNameList[count-1];

    var IA_Path = IA_TRFList;

    var FA_Path = FA_TRFList;



    // Create VTC from FMR in TAL Space

    docVMR.CreateVTCInACPCSpace(InFmrPath,

                                IA_Path,  // IA.trf

                                FA_Path,  // FA.trf

                                ACPC_TRFList,  // *_ACPC.trf

                                OutVTCName,  // output name

                                2,   // datatype integer 2-byte format:1 or float format:2

                                2,   // target resolution 1x1x1:1 or 2x2x2:2 or 3x3x3:3

                                2,   // nearest neighbour:0 or trilinear:1 or sinc:2

                                100  // threshold(Default value:100)

                                );


docVMR.Close();
    bvqx.PrintToLog(OutVTCName + " is done.");

};



bvqx.PrintToLog("Finished.");
