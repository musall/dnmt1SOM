# DNMT1 SOM data analaysis
This repo contains a Matlab-based code collection for  analysis of electrophysiological data in the study by Reichard et al, 2024.
DNMT1-Mediated Regulation of Inhibitory Interneuron Migration Impacts Cortical Architecture and function
https://doi.org/10.1101/2024.09.04.611268

The data for the analysis can be found here: https://doi.org/10.6084/m9.figshare.28282838.v2

The dataset contains Neuropixels recordings from 4 mice that are either SOM-Cre or SOM-DNMT1-KO animals. Each recording contains pre-processed LFP and spiking data that was used to create the electrophysiological results of the study.

To get to the results in the paper, download the dataset with the data in the folder 'Ephys_Data'.
You can then use the main functions (all starting with ds_) to create different figures. For the code to work, the variable 'opts.localPath' needs to be changed to the download Ephys_data folder.

Here is a short description of the included functions:

ds_checkEpileptiformEvents - Test the LFP data for the occurence of epileptiform events and plot their rate and duration

ds_checkOscillatoryPower - Show results for analysis of LFP data as shown in Fig. 7 of the manuscript

ds_showTactileResponse - Show current source density analysis in S1 recordings in response to tactile stimulation

ds_sensoryStimulation - Show results of spiking analysis in response to tactile stimulation

ds_somStimulation - Show neural response to optogenetic stimulation of SOM neurons in S1 and V1
