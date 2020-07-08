READ ME 
by: T.T. Fernandes, September 2019
LarSys - Instituto Superior Tecnico - Universidade de Lisboa

Description

This is a Script that predicts the signal-to-noise ratio (SNR) observed for
different brain tissues (gray matter, white matter, cerebrospinal fluid - CSF) and the contrast-to-noise ratio (CNR) of an acute stroke
lesion relative to those tissues. The implemented code predicts the 
SNR and CNR per time unit for a spin-echo diffusion-weighted sequence, 
taking into account the steady-state value for the longitudinal
magnetization. For that purpose, we adapted an expression used to
predict the SNR per tissue per time unit for a spoiled gradient echo
sequence (Marques, JP et al, 2019). This simulation considers the possibility to use either EPI or spiral readouts and considers the impact of: B0, max gradient amplitude, spatial resolution and b-value (comparing the achieved SNR with that of a typical 1.5T clinical scanner).

How to run the SNR & CNR Study Simulation?

1 - Download the folder of 'jMRI_Publication_Functions'
2 - Select a directory and create a folder called 'SNR_CNR_Study'
3 - Copy files from 'jMRI_Publication_Functions' to this new folder
4 - Run file 'runTests_jMRI_Publication_Functions.m' to perform unit test on all function and scripts of this toolbox
5 - Open file 'SNR_CNR_Study_Simulation.m'
6 - Before running on Matlab, select the values you are interested in see on the outcome plots
7 - Run the script 'SNR_CNR_Study_Simulation.m' - all other functions are called with this script.




