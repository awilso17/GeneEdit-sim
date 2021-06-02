#########################################################################################################################
# Gene flow simulation model for the paper "Gene editing in Farm Animals: A Step Change for Eliminating Epidemics on our Doorstep?" by Petersen et al. 
# Pig breeding pyramid stochastic simulation model - script to run simulation replicates
# calls functions from "GeneFlowModelFcts.R"
# Epidemiological threshold parameters are defined based on outputs from EpidemiologicalModel.R
# AbacusBio Ltd 2019 
#
########################################################################################################################

#!/bin/env Rscript
#setwd("D:/5441")
#setwd("C:/Users/andre/CDocuments/PostDocs/Jaap/Epimodel_paper/PNAS submission")
source("GeneFlowModelFcts.R")

# libraries
library(plyr)
library(dplyr)
library(purrr)

######## Population parameters  ##############################################################################################

indexSD <- 10  # Standard deviation of the economic merit index, using an SD of 1 allows us to scale, but for the results in the report
# a figure of 10 euro was applied. Older Canadian estimates available are
# Yorkshire dam line= 17.2, Landrace = 11.8, Duroc = 9.6 https://www.ccsi.ca/Genetics/Questions%20and%20Answers%20.pdf 
# PropGilts <- 0.242    # proportion of gilts in the sow herd - this is taken care of in the age distribution 

# number of herds per tier
nSPFNuc_Herds <- 3
nProdNuc_Herds <- 10
nMult_Herds <- 20
nBW_Herds <- 300

# age = 4 week batches = months
F_ageFirstMate <- 8
M_ageFirstMate <- 8
# Number of batches between use for mating in sows
FarrowInt <- 5
# remainder for Modulo operator to get around gestation (no double-pregnancy or culling of pregnant sows)
rem <- F_ageFirstMate - FarrowInt
# LitterSize, currently one fixed parameter - could be different for the different SPF breeds, or varied by tier 
LitterSize <- 12

# Selection Props applied (These can be updated to from the gilt requirements model)
# First tier (SPF nucleus)
SPFNucMales_SPFNucFem <- 0.02 # Proportion of male SPF nucleus male offspring retained per batch to mate with SPF females
SPFGiltRepls <- 0.1 # Proportion of female SPF nucleus female offspring retained per batch to mate with SPF males

# Second tier (Production nucleus)
SPFNucMales_ProdNucFem <- 0.1 # Proportion of SPF mat breed A nucleus male offspring per batch transferred to the Production nucleus 
SPFNucFem_ProdNuc <- 0.4 # Proportion of female SPF nucleus female offspring per batch transferred as gilts to the production nucleus
ProdNucGilts_Repls <- 0.2 # Proportion of female production nucleus female offspring retained per batch to enter sow herd

# Third tier (Multiplier)
SPFMales_Mult <- 0.1 # Proportion of SPF mat breed B nucleus male offspring per batch transferred to the multiplier tier 
ProdNucFem_Mult<- 0.5 # Proportion of female production nucleus female offspring per batch transferred as gilts to the multiplier tier

# Fourth tier (Breeder Weaner)
F1Gilts_BW <- 0.6 # Proportion of F1 female offspring from the production nucleus that are transferred as gilts to the breeder weaner tier
TermMales_BW <- 0.2 # Proportion of SPF mat breed B nucleus male offspring per batch transferred to the breeder weaner tier 

# distribution of herd types Mat breed A, Mat breed B, Terminal
SPF_propM0 <- 0.4 # Proportion of males generated in the first generation of the base population
DistHerdTypes<-c(0.4,0.3,0.3) # In the generation of the base population, if there are > 3 SPF herds, assign them to Mat breed A

AgeDist <- c(0.242,0.189,0.147,0.115,0.09,0.07,0.055,0.092) # Distribution of sows by parity (i.e. 24% of sow herd are gilts, 19% parity 1, etc )
AgeDist2 <- rep(AgeDist / 5, each = 5) # convert to 40 batches, each age dist is divided by 5 and replicated 5 times to allow for farrowing interval
AgeDist <- cbind(c(1:8), c(0.242, 0.189, 0.147, 0.115, 0.09, 0.07, 0.055, 0.092)) # create AgeDist dataframe with parity for culling function
colnames(AgeDist) <- c("gest", "dist")
# boar age distribution (assumes that boars are, at the longest, used for 4 years (before they get too heavy to cover / jump))
AgeDist_M <- cbind(c(1:48), c(0.104, 0.068, 0.055, 0.051, 0.046, 0.042, 0.039, 0.037, 0.036, 0.034, 0.033, 0.031, 0.029, 0.027, 0.026, 0.024, 0.023, 0.021, 0.02, 0.019, 0.018, 0.017, 0.016, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011, 0.011, 0.01, 0.01, 0.009, 0.009, 0.008, 0.007, 0.006, 0.006, 0.005, 0.005, 0.004, 0.003, 0.003, 0.002, 0.002, 0.001, 0.001, NA))
colnames(AgeDist_M) <- c("age", "dist")


################## Gene editing parameters ################################################################################################
edit_type <- "C" # C =CRISPR (default), Z =ZFN, T =TALEN, P =Perfect/no errors
edit_mode <- "D" # D = deletion/knock-out, A = insertion/addition 
#edit_prop <- 0.01 # based on TGM (default 0.01 = 1%)
edit_trials <- 1 # number of attempts to edit an embryo (-1 = repeat until success)
embryo_trials <- 1 # number of attempts to transfer an edited embryo (-1 = repeat until success)
# edit_sex <- "M" # sex of animals to be edited

# Gene editing failure rate
# fail rate of edit - affects proportion of successful mutants
fail_rate <- ifelse(edit_type == "C", ifelse(edit_mode == "D", 0.19, 0.37),
                    ifelse(edit_type == "Z", ifelse(edit_mode == "D", 0.45, 0.89),
                           ifelse(edit_type == "T", ifelse(edit_mode == "D", 0.40, 0.79), 0.0)))

death_rate <- ifelse(edit_type == "C", ifelse(edit_mode == "D", 0.39, 0.79),
                     ifelse(edit_type == "Z", ifelse(edit_mode == "D", 0.46, 0.92),
                            ifelse(edit_type == "T", ifelse(edit_mode == "D", 0.44, 0.88), 0.0)))




#SPFeditProp <- 0.2
genoT <- c("SPF","PN")
nMonthsForward = 120


# running replicates:

for (r in 1:10){
  # Running the replicates - if we want to compare scenarios with editing different proportions, need to run a burn in population once for each replicate and 
  # then run the different forward simulation options (including a control where prop edits = 0) from the same burn in, then the burn in is reset between 
  # replicates
  print(paste("Replicate ",r))
  burnInPop <- burnIn(indexSD,nSPFNuc_Herds,nProdNuc_Herds,nMult_Herds,nBW_Herds,AgeDist,AgeDist_M,SPFNucMales_SPFNucFem,SPFGiltRepls,SPFNucMales_ProdNucFem,
                      SPFNucFem_ProdNuc,SPFMales_Mult,ProdNucFem_Mult,ProdNucGilts_Repls,F1Gilts_BW,TermMales_BW,F_ageFirstMate,M_ageFirstMate,FarrowInt,LitterSize)
  
  SPFpop <- burnInPop$SPF
  ProdNucPop <- burnInPop$ProdNuc
  MultPop <- burnInPop$Mult
  BWPop <- burnInPop$BWeaner

  ResEditp2 <- ForwardSim(nMonthsForward,SPFpop,ProdNucPop,MultPop,BWPop,burnInPop$SelProps,F_ageFirstMate,M_ageFirstMate,FarrowInt,LitterSize,AgeDist,AgeDist_M,SPFeditProp=0.2,genoT=c("SPF","PN"))
  ResEditp1 <- ForwardSim(nMonthsForward,SPFpop,ProdNucPop,MultPop,BWPop,burnInPop$SelProps,F_ageFirstMate,M_ageFirstMate,FarrowInt,LitterSize,AgeDist,AgeDist_M,SPFeditProp=0.1,genoT=c("SPF","PN"))
  ResEditp05 <- ForwardSim(nMonthsForward,SPFpop,ProdNucPop,MultPop,BWPop,burnInPop$SelProps,F_ageFirstMate,M_ageFirstMate,FarrowInt,LitterSize,AgeDist,AgeDist_M,SPFeditProp=0.05,genoT=c("SPF","PN"))
  ResEdit0 <- ForwardSim(nMonthsForward,SPFpop,ProdNucPop,MultPop,BWPop,burnInPop$SelProps,F_ageFirstMate,M_ageFirstMate,FarrowInt,LitterSize,AgeDist,AgeDist_M,SPFeditProp=0,genoT=c("SPF","PN"))
  
  ResEditp2All <- cbind(seq(nMonthsForward),ResEditp2$pmerit,ResEditp2$eNumbers,ResEditp2$eSuccess,ResEditp2$nRes)
  ResEditp2All$editProp <- 0.2
  ResEditp2All$rep <- r
  
  ResEditp1All <- cbind(seq(nMonthsForward),ResEditp1$pmerit,ResEditp1$eNumbers,ResEditp1$eSuccess,ResEditp1$nRes)
  ResEditp1All$editProp <- 0.1
  ResEditp1All$rep <- r
  
  ResEditp05All <- cbind(seq(nMonthsForward),ResEditp05$pmerit,ResEditp05$eNumbers,ResEditp05$eSuccess,ResEditp05$nRes)
  ResEditp05All$editProp <- 0.05
  ResEditp05All$rep <- r

  ResEditp0All <- cbind(seq(nMonthsForward),ResEdit0$pmerit,ResEdit0$eNumbers,ResEdit0$eSuccess,ResEdit0$nRes)
  ResEditp0All$editProp <- 0
  ResEditp0All$rep <- r 
  
 
  if(r==1){
    # Collate all the results together into a long table
    AllResN <- rbind(ResEditp2All,ResEditp1All,ResEditp05All,ResEditp0All)
    # AllResN <- rbind(ResEditp1All,ResEditp05All,ResEditp0All)
  } else {
    AllResN <- rbind(AllResN,ResEditp2All,ResEditp1All,ResEditp05All,ResEditp0All)
  }
  
}

# add an ID of replicate_generation - this can be used later to merge the control genetic trends with the different editing scenarios
AllResN$id <- paste(AllResN$rep,AllResN[,1],sep="_")
# set the column numbers 
colnames(AllResN) <- c( "id","gen", "SPF_A","SPF_B","SPF_T","PN","Mult","BW","nCommBorn","ED_SPF_A",
                        "ED_SPF_B", "ED_SPF_T", "ED_PN","ED_Mult","ED_BW","ES_SPF_A", "ES_SPF_B", "ES_SPF_T", "ES_PN","ES_Mult", 
                        "ES_BW","pp_SPF","pp_PN","pp_Mult","pp_BW","pp_BWborn","SizeComm","prop_pp_SPF","prop_pp_PN" ,"prop_pp_Mult",
                        "prop_pp_BW","Spare1","Spare2","editProp","rep")

write.csv(AllResN,paste0("Collated_replicate_results_",Sys.Date(),".csv")
