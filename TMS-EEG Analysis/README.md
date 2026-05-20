# AARATEP

The Advanced Artifact Removal for Automated TMS-EEG Data Processing (AARATEP) found here is a modified version of the original. It's been 
modified to accomodate for our study's modified electrode positions, as well as the lack of explicit ground electrode recording data (our EEG
system subtracts ground from all other channels at the hardware level)

The original version can be found here: https://github.com/chriscline/AARATEPPipeline

Instructions for setting up the environment for the pipeline can also be found in the original repository

NOTE: For ICALabel to run at a reasonable speed, you need to setup matlab to work with XCode (macOS) compiler
