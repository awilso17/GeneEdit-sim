# GeneEdit-sim
This repository contains the simulation code and results for the paper "Gene editing in Farm Animals: A Step Change for Eliminating Epidemics on our Doorstep" by G. E. L. Petersen, J. Buntjer, F. Hely, T. Byrne, B.Whitelaw & A. Doeschl-Wilson
https://www.biorxiv.org/content/10.1101/2021.04.19.440533v1

Scripts:
EpidemiologicalModel.R - R-script to predict impact of gene editing and vaccination on disease prevalence; produces the results for Figures 1-3 in the paper. 
GeneFlowModelFcts.R - R-script containing the functions for the gene flow simulation model. These are called in GeneFlowRunSimulations.R
RunGeneFlowSimulations.R - R-script running the gene flow simulation model; produces the results for Figure 5 in the paper, and Table S1
this models simulates birth, death, genetic selection, and transition of animals between pyramid tiers, the gene editing process the top tear of the pig production pyramid, the genotyping and selection processes in subsequent tiers, and tracks the gene flow of resistance alleles through a pig production pyramid


Data - Results:
