########################################################################################################################
# Functions for the stochastic simulation model GeneFlowRunSimulations
# to genereate gene edited animals and track the gene flow of resistance alleles through a pig breeding pyramid
#
# AbacusBio Ltd 2018 
#
########################################################################################################################

# libraries
library(plyr)
library(dplyr)
library(purrr)

######## Population parameters  ##############################################################################################
# These are the main functions to run the burn in and forward simulation for the gene-editing pig breeding pyramid simulation

######## Functions for key simulation steps ##############################################################################################


createOffspringPRRS <- function(sires,dams,nPiglets,label,indexSD,genNo){
  # create offspring from a given set of sires and dams, including genotype for PRRS resistance based on sire and dam genotypes
  # npiglets is the litter size, label can be a text string indicating which tier/generation offspring is from
  # initiate the data structure for the resulting piglets
  nR <- length(dams$ID)*nPiglets
  piglets=setNames(data.frame(matrix(ncol = 11, nrow = nR)), c("ID", "gen","sex","herd","merit","age","sire","dam","breed","fate","geno"))
  meritSD <-  rnorm(nR, 0, indexSD)

  nR <- length(dams$ID)*nPiglets
  piglets=setNames(data.frame(matrix(ncol = 11, nrow = nR)), c("ID", "gen","sex","herd","merit","age","sire","dam","breed","fate","geno"))
  meritSD <-  rnorm(nR, 0, indexSD)
  
  start <- Sys.time()
  # generate a random order of sires the length of ndams x npiglets per dam
  sireOrder <- sample(1:length(sires$ID), length(dams$ID), replace = TRUE)
  sires_ordered <- sires[c(sireOrder),]
  
  # split the sire and dam genotypes into the separate alleles
  sg <- sires_ordered$geno
  s1 <- substr(sg, 1, 1)
  s2 <- substr(sg, 2, 2)
  dg <- dams$geno
  d1 <- substr(dg, 1, 1)
  d2 <- substr(dg, 2, 2)
  
  # Sample which allele of the PRRS resistant genotype is inherited from each parent
  sireAllele <- sample(c(1,2),nR,replace = TRUE)
  damAllele <- sample(c(1,2),nR,replace = TRUE)
  
  piglets$ID <- paste(label,seq(1:nR),sep="")
  piglets$gen <- genNo
  piglets$sex <- sample(c("M", "F"), nR,replace = TRUE)
  piglets$herd <- dams$herd
  piglets$sire <- sires_ordered$ID
  piglets$dam <- dams$ID
  piglets$merit <- 0.5*(sires_ordered$merit+dams$merit)+meritSD
  piglets$breed <- ifelse(sires_ordered$breed == dams$breed, paste(dams$breed), paste(sires_ordered$breed, dams$breed, sep = ""))
  piglets$age <- -3 # modified to allow for 4 month gestation length 
  piglets$fate <- 1
  piglets$geno <- paste0(ifelse(sireAllele==1,s1,s2),ifelse(damAllele==1,d1,d2))
  
  end <- Sys.time()
  end-start
  
 
  return(piglets)
}

# all the different inheritance cases for genos, assuming no gene drive/preferential inheritance of PRRS resistant edits etc 
# case_when(
#   sg=="PP" & dg=="PP"  ~"PP",
#   sg=="Pp" & dg=="Pp" ~ sample(c("PP","Pp","pp"),1,prob=c(0.25,0.5,0.25)),
#   sg=="PP" & dg=="Pp" | sg=="Pp" & dg=="PP" ~ sample(c("PP","Pp"),1,prob=c(0.5,0.5)),
#   sg=="pp" & dg=="Pp" | sg=="Pp" & dg=="pp" ~ sample(c("pp","Pp"),1,prob=c(0.5,0.5)),
#   sg=="PP" & dg=="pp" | sg == "pp" & dg=="PP" ~ "Pp",
#   sg=="pp" & dg=="pp" ~ "pp"
# )


basePopPRRS <- function(n,indexSD,nSPFNuc_Herds,SPF_propM0,AgeDist,DistHerdTypes,M_ageFirstMate){
  # Initialise a pool of genetics for the SPF Nucleus herds - these are initialised first and then the lower tiers are built up from here
  # set up breeds - must have a min of 3 SPF nucleus herds, 1 of each breed
  if(nSPFNuc_Herds<=3){
    nSPFNuc_Herds<-3
    BreedVec <- cbind(rep(1:nSPFNuc_Herds),c("A","B","T"))
  } else  {
    # if there are more than 3 herds, randomly assign breed
    BreedVec <- cbind(rep(1:nSPFNuc_Herds),sample(c("A","B","T"), size = nSPFNuc_Herds, replace = TRUE, prob = DistHerdTypes))
  }
  colnames(BreedVec) <- c("herd","breed")
  
  # start building up the base pop of SPF animals - then we can start flowing these down the tiers
  SPFnucleus <- as.data.frame(cbind(rep(1:n),rep(0,n),sample(c("M","F"),size=n,replace=TRUE,prob=c(SPF_propM0,1-SPF_propM0)),sample(rep(1:nSPFNuc_Herds),n,replace=TRUE),rnorm(n, mean = 0, sd = indexSD)))
  colnames(SPFnucleus) <- c("ID","gen","sex","herd","merit")
  # assign parity/age by sampling from weighted dist - need to go back 
  SPFnucleus$age <- sample(rep(1:length(AgeDist)), size = n, replace = TRUE, prob = AgeDist)

  # all males set to age 1:first mate
  SPFnucleus$age <- ifelse(SPFnucleus$sex=="M",sample(rep(1:M_ageFirstMate), size = n,replace=TRUE),SPFnucleus$age)
  # merge in breed by herd
  SPFnucleus <- merge(SPFnucleus,BreedVec,by="herd",all.x=TRUE)
  SPFnucleus$fate <- 1 # fate 0/1 where 0 = dead, 1 = alive'
  # All base pop animals start out non-resistant 
  SPFnucleus$geno <- "PP"
  SPFnucleus$merit <- as.numeric(as.character(SPFnucleus$merit))
  
  return(SPFnucleus)  
}

ageing <-function(popdata, mortProp=0.025){
  # function to add ageing/random mortality into the population, the age column is incremented for all animals, then sows over the age of 40 
  # months are culled, and any boars over 36 months are culled. Default mortality is set at 2.5%
  popdata$age <- as.numeric(as.character(popdata$age))
  popdata$age <- popdata$age + 1
  # Check if there are any animals that are too old and need to die, then cull them
  if(nrow(popdata[popdata$sex == "F" & popdata$age > 40 & !is.na(popdata$age),])>0) {popdata[popdata$sex == "F" & popdata$age > 40,]$fate <- 0}
  if(nrow(popdata[popdata$sex == "M" & popdata$age > 36 & !is.na(popdata$age),])>0) {popdata[popdata$sex == "M" & popdata$age > 36,]$fate <- 0}
  
  # add in some random mortality to cut numbers back a little
  popdata$mort <- runif(length(popdata$fate))
  popdata$fate2 <- ifelse(popdata$mort<=mortProp & popdata$fate==1,0,1)
  popdata$mort <- NULL
  popdata[popdata$fate2 == 0,]$fate <- 0
  popdata$fate2 <- NULL
  preage <- nrow(popdata)
  popdata <- subset(popdata,popdata$fate > 0)
  #print(paste0("age die ",preage - nrow(popdata)))
  return(popdata)
}


SowCull <- function(popdata,AgeDist){
  # Function to cull sows to maintain an age parity distribution - sow cull is run in every batch, but only takes into account the sows
  # that are eligible for mating in that match (i.e. pregnant or suckling sows are not considered)
  # split population: 1. check sows vs gilts+males, 2. only select sows that are not currently pregnant
  popdataSows <- subset(popdata, popdata$sex == "F" & popdata$age >= F_ageFirstMate & popdata$age %% FarrowInt == rem & popdata$fate == 1)
  # check if any of these sows are available for culling
  if (length(popdataSows$ID) > 0) {
    # define other part of the population (if there's no sows available for culling, we don't need this..)
    popdataOth <- subset(popdata, !(popdata$sex == "F" & popdata$age >= F_ageFirstMate & popdata$age %% FarrowInt == rem & popdata$fate == 1))
    nOthSows <- nrow(subset(popdataOth, (popdata$sex == "F" & popdata$age >= F_ageFirstMate)))
    # determine parity, store in temporary column
    popdataSows$gest <- ((popdataSows$age - F_ageFirstMate) / FarrowInt)+1
    ageList <- popdataSows %>% group_by(gest) %>% summarise(cnt=n())
    tot <- sum(ageList[, 2])
    
    # merge in the number that should be kept by each age bracket
    ageList <- merge(ageList,AgeDist,by="gest",all.x=TRUE)
    ageList$keep <- round(ageList$dist*tot,0)
    ageList$cnt <- NULL
    ageList$dist <- NULL
    popdataSows <- merge(popdataSows,ageList,by="gest",all.x=TRUE)
       
    # rank by merit within age group, and set the fate of any unrequired sows to -2 
    popdataSows$rank <- unlist(with(popdataSows,tapply(desc(merit),gest,rank)))
    # call fate "2" in case we want to look at culls
    popdataSows$fate <- ifelse(popdataSows$rank>popdataSows$keep,-2,1)
    popdataSows$rank <- NULL
    popdataSows$gest <- NULL
    popdataSows$keep <- NULL
    
    # reassemble data set and return
    popdata <- rbind(popdataSows,popdataOth)
    precull <- nrow(popdata)
    
    popdata <- subset(popdata,popdata$fate>0)
    #print(paste0("Culled: ",precull-nrow(popdata)))
    return(popdata)
  } else { return(popdata) }
}

BoarCull <- function(popdata, AgeDist_M) {
  # Function to keep boars within initial age parity distribution
  ageList <- popdata[popdata$sex=="M" & popdata$fate==1 ,] %>% group_by(age) %>% summarise(cnt <- n())
  #tot <- pmin(sum(ageList[ageList$age>0,2]),maxBoars)
  tot <- sum(ageList[ageList$age>0,2])
  
  # merge in the number that should be kept for each age bracket
  ageList <- merge(ageList,AgeDist_M,by="age",all.x=TRUE)
  ageList$keep <- round(ageList$dist*tot,0)
  ageList$`cnt <- n()` <- NULL
  ageList$dist <- NULL
  popdata <- merge(popdata,ageList,by="age",all.x=TRUE)
  
  # Divide into males vs gilts/sows  
  popdataBoars <- subset(popdata,popdata$sex == "M" & popdata$fate == 1 & !is.na(popdata$keep))
  popdataOth <- subset(popdata,!(popdata$sex == "M" & popdata$fate == 1 & !is.na(popdata$keep)))
  
  # rank by merit within age group, and set the fate of any unrequired boars to -2 
  popdataBoars$rank <- unlist(with(popdataBoars,tapply(desc(merit),age,rank)))
  # call fate "2" in case we want to look at culls
  popdataBoars$fate <- ifelse(popdataBoars$rank>popdataBoars$keep,-2,1)
  popdataBoars$rank <- NULL
  
  # reassemble data set and return
  popdata <- rbind(popdataBoars,popdataOth)
  popdata$keep <- NULL
  
  popdata <- subset(popdata,popdata$fate>0)
  return(popdata)
}

# Gene editing function
edit_genes <- function(animals, edit_prop = 0.0, edit_type, edit_trials, embryo_trials, fail_rate, death_rate, sel_prop) { 
  # gene editing workflow
  #     0.0. Subset individuals that would have to be edited
  #     0.1. Sort individuals based on merit to identify the top animals to edit
  # For each recessive to be edited:
  #     1. Select the top edit_prop proportion of animals, currently at least 1 animal always will be edited
  #     2. Do the edit for PP, Pp and pP genotypes
  #     3. Check to see if the edit succeeded
  #     4. Update the animal's genotype
  #     5. Update the edit_status list
  #     6. Update the ge_trials (number of tries until successful edit)
  
  ## 0. Subset individuals that would have to be edited
  # Check and see if we have recessives to edit.
  homozygotes <- subset(animals, animals$geno == "pp")
  others <- subset(animals, animals$geno == "PP" | animals$geno == "pP" | animals$geno == "Pp")
  # check that homozygotes do not already supply all the required animals for selection
  if (length(homozygotes$ID) < ceiling(length(animals$ID) * sel_prop)) {
      n_edit <- ceiling(length(animals$ID) * edit_prop)
      if (n_edit >= 1) {
            cat("Attempting to edit ", n_edit, " animals edited with edit_prop = ", edit_prop, "and edit_type = ", edit_type)
        }
      ## 0.1. Sort individuals based on merit ##
      others <- others[order(others$merit, decreasing = TRUE),]
      # add columns to allow for recording editing information
      others$edit_status <- 0
      others$ge_trials <- others$et_trials <- NA
  
      ## 1. Select the top edit_prop proportion of animals.
      edits <- 0
      for (animal in 1:length(others$ID)) {
        if(edits < n_edit){
          # assuming that both heterozygotes and undesirable homozygotes need to be edited (and desirable homozygotes do not):
          # if (others$geno[animal] == "PP" | others$geno[animal] == "Pp" | animals$geno[animal] == "pP") {  - Don't need this anymore, as homozygotes are already removed
            # 2. Do the edit for PP and Pp genotypes
            if (edit_trials > 0) {
              edits <- edits+1
              # 3a. (i) if edit_trials > 0 then only a fixed number of trials will be carried out. 
              # If there is no success before the final trial then the editing process fails.
              outcomes <- rbinom(edit_trials, 1, (1 - fail_rate))
              # check if edit was successful
              if (1 %in% outcomes) {
                # 4. Update the animal's genotype - for now assume that we can edit both in the same editing step
                #animals$geno[animal] <- ifelse((animals$geno[animal] == "Pp"), "PP", "pp")
                others$geno[animal] <-  "pp"
                # 5. Update the edit_status 
                others$edit_status[animal] <- 1
                # 6. Update the animal's edit count with the time of the first successful edit
                others$ge_trials[animal] <- which(outcomes != 0)[1]
            
                # 3b. Was the successfully edited embryo carried to term?
                if (embryo_trials > 0){
                  # 3b. (i) If the embryo died then we need to update the fate. If edit_trials > 0 then only a fixed number of trials
                  #         will be carried out. If there is no success before the final trial then the editing process
                  #         fails.
                  outcomes <- rbinom(edit_trials, 1, (1 - death_rate))
                  if (any(outcomes)){
                    others$et_trials[animal] <- which(outcomes != 0)[1]
                  } else {
                    # If ET fails, animal is dead/not available for selection
                    others$et_trials[animal] <- embryo_trials
                    others$fate[animal] <- 0
                  } 
                }
              } else {
                others$edit_status[animal] <- -1
              }
            }
            else { print ("edit_trials should never be 0, skipping editing step!") }
          #}
        }
      }
      # print how many animals were edited AND survived the ET
      #print(paste(sum(others$edit_status == 1 & others$fate == 1), "animals successfully edited", sep = " "))
      # record the number of edits attempted, and the number of successful edits to an output list
      edits <- sum(others$edit_status == 1) + sum(others$edit_status == -1)
      success <- sum(others$edit_status == 1 & others$fate == 1)
      animals <- rbind.fill(homozygotes, others)
      return(list(animals=animals, edits = edits, success = success))
  }
  else { cat("Enough resistant animals in the population for selection, no individuals edited")
    # If we have enough homozygotes available for to meet the selection 
      edits <- 0
      success <- 0
      animals$edit_status <- NA
      animals$et_trials <- NA
      animals$ge_trials <- NA
      animals <- animals
      return(list(animals = animals, edits = edits, success = success))
  }
}


# Selecting top gilts/boars to retain based on homozygote PRRS resistant first (similar to parent selection without age requirements)
retainGeno <- function(pop, nsel) {
    # rearrange population based on "increasing" genotype (lower case comes before UPPER CASE) and decreasing merit
    pop <- subset(pop, pop$fate == 1)
    selList <- pop[order(pop$geno, - as.numeric(pop$merit)),]
    # Take the top nSel in order
    selList <- selList[1:nsel,]
    return(selList)
}

# Run a burn in to set up the rest of the pyramid until we have 8 parities of sows in the BW tier
burnIn <- function(indexSD,nSPFNuc_Herds,nProdNuc_Herds,nMult_Herds,nBW_Herds,AgeDist,AgeDist_M,SPFNucMales_SPFNucFem,SPFGiltRepls,SPFNucMales_ProdNucFem,
                   SPFNucFem_ProdNuc,SPFMales_Mult,ProdNucFem_Mult,ProdNucGilts_Repls,F1Gilts_BW,TermMales_BW,F_ageFirstMate,M_ageFirstMate,FarrowInt,LitterSize){
  AgeDistt <- cbind(rep(1:length(AgeDist[,2])),AgeDist[,2])
  colnames(AgeDistt)<-c("age","dist")
  AgeDist2 <-  rep(AgeDist[,2] / 5, each = 5)
  
  MaxSowsSPF <- 6000
  MaxSowsPN <- 33000
  MaxSowsM <- 50000
  MaxSowsBW <- 500000
  
  # create a base population of SPF animals - initial pop of 6000 (2000/breed) is about in line with Visscher's "modern" pyramid with 2000 sows in top tier
  SPFpop <- basePopPRRS(6000,1,nSPFNuc_Herds,SPF_propM0,AgeDist2,DistHerdTypes,M_ageFirstMate)
  SPFpop$merit <- as.numeric(as.character(SPFpop$merit))
  SPFpop$ID <- as.character(SPFpop$ID)
  SPFpop$breed <- as.character(SPFpop$breed)
  
  # When selecting dams, need them every FarrowInt batches after F_ageFirstMate, so age modulo farrowInt == diff between can get around this
  rem <- F_ageFirstMate-FarrowInt 
  
  # Set up production nucleus, multiplier and breeder weaner data frames
  ProdNucPop <- data.frame(matrix(ncol=length(names(SPFpop)),nrow=0))
  colnames(ProdNucPop) <- names(SPFpop)
  
  MultPop <- data.frame(matrix(ncol=length(names(SPFpop)),nrow=0))
  colnames(MultPop) <- names(SPFpop)
  
  BWPop <- data.frame(matrix(ncol=length(names(SPFpop)),nrow=0))
  colnames(BWPop) <- names(SPFpop)
  
  # run for 70 "burn in" batches to set up production nucleus
  for(g in -70:0){
    print(paste("Gen:",g,sep=" "))
    
    # Maternal breed A
    # males must be "age at first mating" or older, and have fate = 1
    matAmales <- subset(SPFpop,SPFpop$sex=="M" & SPFpop$breed=="A" & SPFpop$age >= M_ageFirstMate & SPFpop$fate==1)
    # dams from breed A, greater than age of first mating, and mod farrow int == 3 (a way of avoiding pregnant sows)
    matAdams <- subset(SPFpop,SPFpop$sex=="F" & SPFpop$breed=="A" & SPFpop$age>= F_ageFirstMate & SPFpop$age %% FarrowInt==rem & SPFpop$fate==1)
    
    # Create new gen of maternal breed A
    newGenASPF <- createOffspringPRRS(matAmales,matAdams,LitterSize,paste0("SPF_",g,"_"),indexSD,g)
    
    # Subset males and females out and sort by decreasing merit 
    newGenASPF_F <- subset(newGenASPF,newGenASPF$sex=="F")
    newGenASPF_F <- newGenASPF_F[order(newGenASPF_F$merit,decreasing=TRUE),]
    newGenASPF_M <- subset(newGenASPF,newGenASPF$sex=="M")
    newGenASPF_M <- newGenASPF_M[order(newGenASPF_M$merit,decreasing=TRUE),]
    
    # Repeat for maternal breed B of SPF tier - select eligible parents, create a new generation
    matBmales <- subset(SPFpop,SPFpop$sex=="M" & SPFpop$breed=="B" & SPFpop$age >= M_ageFirstMate & SPFpop$fate==1)
    matBdams <- subset(SPFpop,SPFpop$sex=="F" & SPFpop$breed=="B" & SPFpop$age>= F_ageFirstMate & SPFpop$age %% FarrowInt==rem & SPFpop$fate==1)
    
    newGenBSPF <- createOffspringPRRS(matBmales,matBdams,LitterSize,paste0("SPF_",g,"_"),indexSD,g)
    
    newGenBSPF_F <- subset(newGenBSPF,newGenBSPF$sex=="F")
    newGenBSPF_F <- newGenBSPF_F[order(newGenBSPF_F$merit,decreasing=TRUE),]
    
    newGenBSPF_M <- subset(newGenBSPF,newGenBSPF$sex=="M")
    newGenBSPF_M <- newGenBSPF_M[order(newGenBSPF_M$merit,decreasing=TRUE),]
    
    # Repeat for terminal breed T of SPF tier - select eligible parents, create a new generation
    matTmales <- subset(SPFpop,SPFpop$sex=="M" & SPFpop$breed=="T" & SPFpop$age >= M_ageFirstMate & SPFpop$fate==1)
    matTdams <- subset(SPFpop,SPFpop$sex=="F" & SPFpop$breed=="T" & SPFpop$age>= F_ageFirstMate & SPFpop$age %% FarrowInt==rem & SPFpop$fate==1)
    
    newGenTSPF <- createOffspringPRRS(matTmales,matTdams,LitterSize,paste0("SPF_",g,"_"),indexSD,g)
    
    newGenTSPF_F <- subset(newGenTSPF,newGenTSPF$sex=="F")
    newGenTSPF_F <- newGenTSPF_F[order(newGenTSPF_F$merit,decreasing=TRUE),]
    newGenTSPF_M <- subset(newGenTSPF,newGenTSPF$sex=="M")
    newGenTSPF_M <- newGenTSPF_M[order(newGenTSPF_M$merit,decreasing=TRUE),]
    
    if(g==-70){
      # In first gen of the burn in establish the SPF selection numbers = selection proportion x number of gilts/boars available
      selnF_SPF <- round(SPFGiltRepls*length(newGenASPF_F[,1]))
      selnM_SPF <- round(SPFNucMales_SPFNucFem*length(newGenASPF_M[,1]))
      selnF_PN <- round(SPFNucFem_ProdNuc*length(newGenASPF_F[,1]))
      selnM_PN <- round(SPFNucMales_ProdNucFem*length(newGenASPF_M[,1]))
      #print(paste("SPF",selnF_SPF,selnM_SPF,selnF_PN,selnM_PN,sep=","))
    } else {
      # to avoid the size of the SPF herd ballooning, for each subsequent generation we take the minimum of 
      # selectionProp*nOffspring and the previous number selected
      selnF_SPF <- min(selnF_SPF,round(SPFGiltRepls*length(newGenASPF_F[,1])))
      selnM_SPF <- min(selnM_SPF,round(SPFNucMales_SPFNucFem*length(newGenASPF_M[,1])))
      selnF_PN <- min(selnF_PN,round(SPFNucFem_ProdNuc*length(newGenASPF_F[,1])))
      selnM_PN <- min(selnM_PN,round(SPFNucMales_ProdNucFem*length(newGenASPF_M[,1])))
    }
    
    # start a list for all the selection proportions to return at the end of the loop - this is then be used as an input to the 
    # forward simulation
    selNumbers <- list(SPFgilts=selnF_SPF,SPFboars=selnM_SPF,SPFpngilts=selnF_PN,SPFpnboars=selnM_PN)
    
    # Add the n selected gilts and boars from offspring to the SPF for each of the three breeds using the selection numbers defined above
    SPFpop <- rbind.fill(SPFpop,newGenASPF_F[1:selnF_SPF,],newGenASPF_M[1:selnM_SPF,],
                         newGenBSPF_M[1:selnM_SPF,],newGenBSPF_F[1:selnF_SPF,],newGenTSPF_M[1:selnM_SPF,],newGenTSPF_F[1:selnF_SPF,])
    
    #print the number of live males/females in the SPF pop data set 
    print(paste0("Pop stats: M ",sum(SPFpop$sex=="M" & SPFpop$fate==1),", F ",sum(SPFpop$sex=="F" & SPFpop$fate==1)))
    
    # Add the next best n selected gilts/boars from the SPF mat breed A into the production nucleus 
    ProdNucPop <-rbind.fill(ProdNucPop,newGenASPF_F[(selnF_SPF+1):(selnF_SPF+selnF_PN),],
                       newGenASPF_M[(selnM_SPF+1):(selnM_SPF+selnM_PN),])
    
    # Once we have 30 months (3 age classes of sows) in the Prod nucleus, start rolling them forward and building Mult tier too  
    if (g >= -40){
      # Some housekeeping as merit/age kept getting converted to a character and this makes sure they stay numeric
      ProdNucPop$merit <- as.numeric(as.character(ProdNucPop$merit))
      ProdNucPop$ID <- as.character(ProdNucPop$ID)
      ProdNucPop$breed <- as.character(ProdNucPop$breed)
      ProdNucPop$age <- as.numeric(as.character(ProdNucPop$age))
      
      # Subset the available boars/sows from in the production nucleus
      PN_males <- subset(ProdNucPop,ProdNucPop$sex=="M" & ProdNucPop$breed=="A" & ProdNucPop$age >= M_ageFirstMate & ProdNucPop$fate==1)
      PN_dams <- subset(ProdNucPop,ProdNucPop$sex=="F" & ProdNucPop$breed=="A" & ProdNucPop$age >= F_ageFirstMate & ProdNucPop$age %% FarrowInt==rem & ProdNucPop$fate==1)
      
      # create new generation of prod nucleus, top gilts will be retained, and gilts transferred to the multiplier
      newGenPN <- createOffspringPRRS(PN_males,PN_dams,LitterSize,paste0("PN_",g,"_"),indexSD,g)
      newGenPN_F <- subset(newGenPN,newGenPN$sex=="F")
      newGenPN_F <- newGenPN_F[order(newGenPN_F$merit,decreasing = TRUE),]
      
      if(g==-40){
        # establish selection numbers
        selnF_PNG <- round(ProdNucGilts_Repls*length(newGenPN_F[,1]))
        selnF_PN_M <- round(ProdNucFem_Mult*length(newGenPN_F[,1]))
        selnM_B_M <- round(SPFMales_Mult*length(newGenBSPF_M[,1]))
        #print(paste("PN",selnF_PNG,selnF_PN_M,selnM_B_M,sep=","))
      } else {
        selnF_PNG <- min(selnF_PNG,round(ProdNucGilts_Repls*length(newGenPN_F[,1])))
        selnF_PN_M <- min(selnF_PN_M,round(ProdNucFem_Mult*length(newGenPN_F[,1])))
        selnM_B_M <- min(selnM_B_M,round(SPFMales_Mult*length(newGenBSPF_M[,1])))
      }
      # Add these selection numbers to the list 
      selNumbers <- c(selNumbers,list(PNgilts=selnF_PNG,PNmgilts=selnF_PN_M,SPFmboars=selnM_B_M))
      
      # Add the retained gilts from the production nucleus to the PN pop
      ProdNucPop <-rbind.fill(ProdNucPop,newGenPN_F[1:selnF_PNG,])
      
      # every FarrowInt generations run the sow cull 
      if(g>= -25){
        # let 15 batches of gilts be transferred before we start trying to cull 
        ProdNucPop <- SowCull(ProdNucPop,AgeDist)
      }
      
      # Start funneling the production nucleus gilts and Maternal breed B males into the multiplier tier
      MultPop <- rbind.fill(MultPop,newGenPN_F[(selnF_PNG+1):(selnF_PNG+selnF_PN_M),],newGenBSPF_M[(selnM_SPF+1):(selnM_SPF+selnM_B_M),])
      
      # run the age function for the multiplier population 
      MultPop <- ageing(MultPop)
      
      if(g>=-20){
        # Once there are 50 batches of SPF, Prod Nuc and Mult, start creating the breeder weaner tier
        MultPop$merit <- as.numeric(as.character(MultPop$merit))
        MultPop$ID <- as.character(MultPop$ID)
        MultPop$breed <- as.character(MultPop$breed)
        MultPop$age <- as.numeric(as.character(MultPop$age))
        
        # Select available sows and boars from the multiplier tier 
        Mult_males <- subset(MultPop,MultPop$sex=="M" & MultPop$breed=="B" & MultPop$age >= M_ageFirstMate & MultPop$fate == 1)
        Mult_dams <- subset(MultPop,MultPop$sex=="F" & MultPop$age >= F_ageFirstMate & MultPop$age %% FarrowInt==rem & MultPop$fate==1)
        
        # create new generation of prod nucleus, top gilts will be retained, and gilts transferred to the multiplier
        newGenMult <- createOffspringPRRS(Mult_males,Mult_dams,LitterSize,paste0("M_",g,"_"),indexSD,g)
        newGenMult_F <- subset(newGenMult,newGenMult$sex=="F")
        newGenMult_F <- newGenMult_F[order(newGenMult_F$merit,decreasing = TRUE),]
        
        if(g==-20){
          # establish selection numbers
          selnM_T_BW <- round(TermMales_BW*length(newGenTSPF_M[,1]))
          selnF_M_BW <- round(F1Gilts_BW*length(newGenMult_F[,1]))
          #print(paste("BW:",selnF_M_BW,selnM_T_BW,sep=","))
        } else {
          selnM_T_BW <- min(selnM_T_BW,round(TermMales_BW*length(newGenTSPF_M[,1])))
          selnF_M_BW <- min(selnF_M_BW,round(F1Gilts_BW*length(newGenMult_F[,1])))
        }
        # add the selection numbers for the BW tier to the list
        selNumbers <- c(selNumbers,list(Mgilts=selnF_M_BW,SPFtboars=selnM_T_BW))
        
        
        #BWPop <- rbind.fill(BWPop,newGenMult_F[1:selnF_M_BW,],newGenTSPF_M[(selnM_SPF+1):(selnM_SPF+selnM_T_BW),])
        # In the burn in phase to get the numbers up faster we take all available gilts from the multiplier to establish the BW sow population  
        BWPop <- rbind.fill(BWPop,newGenMult_F,newGenTSPF_M[(selnM_SPF+1):(selnM_SPF+selnM_T_BW),])
        BWPop <- ageing(BWPop)
        
        BWPop$age <- as.numeric(as.character(BWPop$age))
        BWPop$merit <- as.numeric(as.character(BWPop$merit))
        BWPop$breed <- as.character(BWPop$breed)
        
        
        # only start sow cull in Multiplier tier once we are mating within the tier
        if(g >= -5){
          MultPop <- SowCull(MultPop,AgeDist)
        }
      }
        
      # Add the retained gilts from the production nucleus to the PN pop
      ProdNucPop <-rbind.fill(ProdNucPop,newGenPN_F[1:selnF_PNG,])
      
      # Start funneling the ProdNuc gilts and Maternal breed B males into the multiplier tier
      MultPop <- rbind.fill(MultPop,newGenPN_F[(selnF_PNG+1):(selnF_PNG+selnF_PN_M),],newGenBSPF_M[(selnM_SPF+1):(selnM_SPF+selnM_B_M),])
      
      
    }
    

    # Run the sow cull function for the SPF
    SPFpop <- SowCull(SPFpop,AgeDist)
    
    # Run the ageing function for the SPF
    SPFpop <- ageing(SPFpop)
    SPFpop$merit <- as.numeric(as.character(SPFpop$merit))
    SPFpop$ID <- as.character(SPFpop$ID)
    
    # Run the ageing function for the production nucleus
    ProdNucPop <- ageing(ProdNucPop)
 
  }

  # add population name identifier to each data frame to show which population they are in 
  SPFpop$popname <- "SPFNUC"
  ProdNucPop$popname <- "PRODNUC"
  MultPop$popname <- "MULT"
  BWPop$popname <- "BW"
  
  #Combine tiers for output 
  AllPop <- rbind.fill(SPFpop,ProdNucPop,MultPop,BWPop)
  
  # If required, write the burn in population to file for faster testing
  # write.csv(AllPop,paste0("TestPopulationBurnin_functest",Sys.Date(),".csv"),row.names = FALSE,quote = FALSE)
  
  return(list(SPF=SPFpop,ProdNuc=ProdNucPop,Mult=MultPop,BWeaner=BWPop,SelProps=selNumbers)) 
}

# testing reset - if you want to reset to the burnin go back to here #
# SPFpop <- subset(AllPop,popname=="SPFNUC")
# SPFpop$merit <- as.numeric(as.character(SPFpop$merit))
# SPFpop$age <- as.numeric(as.character(SPFpop$age))
# SPFpop$ID <- as.character(SPFpop$ID)
# ProdNucPop <- subset(AllPop,popname=="PRODNUC")
# MultPop<- subset(AllPop,popname=="MULT")
# BWPop<- subset(AllPop,popname=="BW")

############################################################################################################################################
#Start of the rolling forward simulation
############################################################################################################################################

# This function allows replicates to be run fromt he same burn in population - i.e. comparing the effect of editing different proportions
ForwardSim <- function(timeSteps,SPFpop,ProdNucPop,MultPop,BWPop,SelNumbers,F_ageFirstMate,M_ageFirstMate,FarrowInt,LitterSize,AgeDist,AgeDist_M,SPFeditProp=0.1,genoT=c("SPF","PN")){

  # Unpack the selection numbers from the burn in output list
  selnF_SPF <- SelNumbers$SPFgilts 
  selnM_SPF <- SelNumbers$SPFboars
  selnF_PN <- SelNumbers$SPFpngilts
  selnM_PN <- SelNumbers$SPFpnboars
  selnM_B_M <- SelNumbers$SPFmboars
  selnM_T_BW <- SelNumbers$SPFtboars
  selnF_PNG <- SelNumbers$PNgilts
  selnF_PN_M <- SelNumbers$PNmgilts
  selnF_M_BW <- SelNumbers$Mgilts
  
  # initialise data frames for the commerical piglets, along with some data frames to hold results 
  CommPop <- data.frame(matrix(ncol=length(names(SPFpop)),nrow=0))
  
  pigletmerit <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(pigletmerit) <- c("SPF_A", "SPF_B", "SPF_T", "PN", "Mult", "BW","commBrn")
  edit_numbers <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(edit_numbers) <- c("SPF_A", "SPF_B", "SPF_T", "PN", "Mult", "BW")
  edit_success <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(edit_success) <- c("SPF_A", "SPF_B", "SPF_T", "PN", "Mult", "BW")
  resistant_vec <- data.frame(matrix(ncol = 6, nrow = 0))
  propres_vec <- data.frame(matrix(ncol = 6, nrow = 0))
  
  # main loop in the forward simulation 
  for(g in rep(1:timeSteps)){
    print(paste("Iteration: ", g))
    print(Sys.time())
    # Initialise vectors to hold the results - these will be appended to the correct data frame
    pmvec <- c(0, 0, 0, 0, 0, 0)
    editvec <- c(0, 0, 0, 0, 0, 0)
    successvec <- c(0, 0, 0, 0, 0, 0)
    numbervec <- c(0, 0, 0, 0, 0, 0)
    propresvec <- c(0, 0, 0, 0, 0, 0)


    #### SPF population ####
    # run cull BEFORE animals are selected as potential parents
    SPFpop <- SowCull(SPFpop, AgeDist)
    SPFpop <- BoarCull(SPFpop, AgeDist_M)


    ## Create new generation of Maternal breed A, record mean merit, and gene editing results ##
    matAmales <- subset(SPFpop,SPFpop$sex=="M" & SPFpop$breed=="A" & SPFpop$age >= M_ageFirstMate & SPFpop$fate==1)
    # dams from breed A, greater than age of first mating, and mod farrow int == 3 (a way of avoiding "double mating" of already pregnant sows)
    matAdams <- subset(SPFpop,SPFpop$sex=="F" & SPFpop$breed=="A" & SPFpop$age >= F_ageFirstMate & SPFpop$age %% FarrowInt==rem & SPFpop$fate==1)
    
    newGenASPF <- createOffspringPRRS(matAmales,matAdams,LitterSize,paste0("SPF",g,"_"),indexSD,g)
    # Take the average merit of piglets born in the new generation 
    pmvec[1] <- mean(as.numeric(as.character(newGenASPF$merit)), na.rm = TRUE)

    newGenASPF_F <- subset(newGenASPF,newGenASPF$sex=="F")
    newGenASPF_F <- newGenASPF_F[order(newGenASPF_F$merit,decreasing=TRUE),]

    newGenASPF_M <- subset(newGenASPF, newGenASPF$sex == "M")
    newGenASPF_M <- newGenASPF_M[order(newGenASPF_M$merit, decreasing = TRUE),]

    # editing in maternal breed A - both males and females
    editASPF_F <- edit_genes(newGenASPF_F, SPFeditProp, edit_type, edit_trials, embryo_trials, fail_rate, death_rate, (SPFNucFem_ProdNuc + ProdNucGilts_Repls))
    editASPF_M <- edit_genes(newGenASPF_M, SPFeditProp, edit_type, edit_trials, embryo_trials, fail_rate, death_rate, SPFNucMales_ProdNucFem)
    # Record the total number of attempted edits, and the number that were successful 
    editvec[1] <- editASPF_F$edits+editASPF_M$edits
    successvec[1] <- editASPF_F$success+editASPF_M$success
    
    newGenASPF_F <- editASPF_F$animals
    newGenASPF_M <- editASPF_M$animals
    
    newGenASPF_F <- newGenASPF_F[order(newGenASPF_F$merit,decreasing=TRUE),]
    newGenASPF_M <- newGenASPF_M[order(newGenASPF_M$merit,decreasing=TRUE),]
  
    ## Create new generation of Maternal breed B, record mean merit, and gene editing results ##
    matBmales <- subset(SPFpop,SPFpop$sex=="M" & SPFpop$breed=="B" & SPFpop$age >= M_ageFirstMate & SPFpop$fate==1)
    matBdams <- subset(SPFpop,SPFpop$sex=="F" & SPFpop$breed=="B" & SPFpop$age >= F_ageFirstMate & SPFpop$age %% FarrowInt==rem & SPFpop$fate==1)
    
    newGenBSPF <- createOffspringPRRS(matBmales,matBdams,LitterSize,paste0("SPF",g,"_"),indexSD,g)
    pmvec[2] <- mean(as.numeric(as.character(newGenBSPF$merit)), na.rm = TRUE)
    
    newGenBSPF_F <- subset(newGenBSPF,newGenBSPF$sex=="F")
    newGenBSPF_F <- newGenBSPF_F[order(newGenBSPF_F$merit,decreasing=TRUE),]
    
    newGenBSPF_M <- subset(newGenBSPF,newGenBSPF$sex=="M")
    
    # Editing in maternal breed B, males only
    editBSPF_M <- edit_genes(newGenBSPF_M, SPFeditProp, edit_type, edit_trials, embryo_trials, fail_rate, death_rate, SPFNucMales_ProdNucFem)
    editvec[2] <- editBSPF_M$edits
    successvec[2] <- editBSPF_M$success
    
    newGenBSPF_M <- editBSPF_M$animals
    newGenBSPF_M <- newGenBSPF_M[order(newGenBSPF_M$merit,decreasing=TRUE),]
    
    
    ## Create new generation of Terminal breed T, record mean merit, and gene editing results ##
    matTmales <- subset(SPFpop,SPFpop$sex=="M" & SPFpop$breed=="T" & SPFpop$age >= M_ageFirstMate & SPFpop$fate==1)
    matTdams <- subset(SPFpop,SPFpop$sex=="F" & SPFpop$breed=="T" & SPFpop$age >= F_ageFirstMate & SPFpop$age %% FarrowInt==rem & SPFpop$fate==1)
    
    newGenTSPF <- createOffspringPRRS(matTmales,matTdams,LitterSize,paste0("SPF",g,"_"),indexSD,g)
    pmvec[3] <- mean(as.numeric(as.character(newGenTSPF$merit)), na.rm = TRUE)
    
    newGenTSPF_F <- subset(newGenTSPF,newGenTSPF$sex=="F")
    newGenTSPF_F <- newGenTSPF_F[order(newGenTSPF_F$merit,decreasing=TRUE),]
    
    newGenTSPF_M <- subset(newGenTSPF,newGenTSPF$sex=="M")
    # Editing in terminal breed T, males only
    editTSPF_M <- edit_genes(newGenTSPF_M, SPFeditProp, edit_type, edit_trials, embryo_trials, fail_rate, death_rate, SPFNucMales_ProdNucFem)
    editvec[3] <- editTSPF_M$edits
    successvec[3] <- editTSPF_M$success
    
    newGenTSPF_M <- editTSPF_M$animals
    newGenTSPF_M <- newGenTSPF_M[order(newGenTSPF_M$merit,decreasing=TRUE),]
    
    # If we are editing and genotyping (i.e. we will know the PRRS genotypes of animals)
    if("SPF" %in% genoT){
      # select PRRS resistant geno first for all SPF replacements
      matA_gilts <- retainGeno(newGenASPF_F,selnF_SPF)
      matB_gilts <- retainGeno(newGenBSPF_F,selnF_SPF)
      T_gilts <- retainGeno(newGenTSPF_F,selnF_SPF)
      
      matA_boars <- retainGeno(newGenASPF_M,selnM_SPF)
      matB_boars <- retainGeno(newGenBSPF_M,selnM_SPF)
      T_boars <- retainGeno(newGenTSPF_M,selnM_SPF)
      
      # also assume that SPF gilts and boars in lower tiers would be preferentially selected
      SPFgilts_PN <- retainGeno(newGenASPF_F,(selnF_SPF+selnF_PN))[(selnF_SPF+1):(selnF_SPF+selnF_PN),]
      SPFboars_PN <- retainGeno(newGenASPF_M,(selnM_SPF+selnM_PN))[(selnM_SPF+1):(selnM_SPF+selnM_PN),]
      
      Mult_boars <- retainGeno(newGenBSPF_M,(selnM_SPF+selnM_B_M))[(selnM_SPF+1):(selnM_SPF+selnM_B_M),]
      BW_boars <- retainGeno(newGenTSPF_M,(selnM_SPF+selnM_T_BW))[(selnM_SPF+1):(selnM_SPF+selnM_T_BW),]
      
    } else {
      # Otherwise select top merit animals regardless of genotype
      matA_gilts <- newGenASPF_F[1:selnF_SPF,]
      matB_gilts <- newGenBSPF_F[1:selnF_SPF,]
      T_gilts <- newGenTSPF_F[1:selnF_SPF,]
      
      matA_boars <- newGenASPF_M[1:selnM_SPF,]
      matB_boars <- newGenBSPF_M[1:selnM_SPF,]
      T_boars <- newGenTSPF_M[1:selnM_SPF,]
      
      SPFgilts_PN <- newGenASPF_F[(selnF_SPF+1):(selnF_SPF+selnF_PN),]
      SPFboars_PN <- newGenASPF_M[(selnM_SPF+1):(selnM_SPF+selnM_PN),]
      
      Mult_boars <- newGenBSPF_M[(selnM_SPF+1):(selnM_SPF+selnM_B_M),]
      BW_boars <- newGenTSPF_M[(selnM_SPF+1):(selnM_SPF+selnM_T_BW),]
    }
    
    # add the top males and females from offspring to the SPF
    SPFpop <- rbind.fill(SPFpop,matA_gilts,matA_boars,matB_gilts,matB_boars,T_gilts,T_boars)
    
    SPFpop <- ageing(SPFpop)
    SPFpop$merit <- as.numeric(as.character(SPFpop$merit))
    SPFpop$age <- as.numeric(as.character(SPFpop$age))
    SPFpop$ID <- as.character(SPFpop$ID)
    
    #print(paste0("SPF pop stats: M ",sum(SPFpop$sex=="M" & SPFpop$fate==1),", F ",sum(SPFpop$sex=="F" & SPFpop$fate==1)))

      
    #### Production nucleus ####
    # run cull BEFORE transferred animals are added and potential parents are selected
    ProdNucPop <- SowCull(ProdNucPop, AgeDist)
    ProdNucPop <- BoarCull(ProdNucPop, AgeDist_M)
    
    # add transferred animals from SPF nucleus
    ProdNucPop <-rbind.fill(ProdNucPop,SPFgilts_PN,SPFboars_PN)
    
    ProdNucPop$merit <- as.numeric(as.character(ProdNucPop$merit))
    ProdNucPop$ID <- as.character(ProdNucPop$ID)
    ProdNucPop$breed <- as.character(ProdNucPop$breed)
    ProdNucPop$age <- as.numeric(as.character(ProdNucPop$age))

    PN_males <- subset(ProdNucPop,ProdNucPop$sex=="M" & ProdNucPop$breed=="A" & ProdNucPop$age >= M_ageFirstMate & ProdNucPop$fate==1)
    PN_dams <- subset(ProdNucPop,ProdNucPop$sex=="F" & ProdNucPop$breed=="A" & ProdNucPop$age >= F_ageFirstMate & ProdNucPop$age %% FarrowInt==rem & ProdNucPop$fate==1)
    
    # create new generation of prod nucleus, top gilts will be retained, and gilts transferred to the multiplier
    newGenPN <- createOffspringPRRS(PN_males,PN_dams,LitterSize,paste0("PN",g,"_"),indexSD,g)
    pmvec[4] <- mean(as.numeric(as.character(newGenPN$merit)), na.rm = TRUE)
    
    newGenPN_F <- subset(newGenPN,newGenPN$sex=="F")
    newGenPN_F <- newGenPN_F[order(newGenPN_F$merit,decreasing = TRUE),]
    
    # If we are editing and genotyping (i.e. we will know the PRRS genotypes of animals)
    if("PN" %in% genoT){
      # select PRRS resistant geno first for all SPF replacements
      PN_gilts <- retainGeno(newGenPN_F,selnF_PNG)
      Mult_gilts <- retainGeno(newGenPN_F,(selnF_PNG+selnF_PN_M))[(selnF_PNG+1):(selnF_PNG+selnF_PN_M),]
    } else {
      # Otherwise select on merit alone
      PN_gilts <- newGenPN_F[1:selnF_PNG,]
      Mult_gilts <- newGenPN_F[(selnF_PNG+1):(selnF_PNG+selnF_PN_M),]
    }
    
    # Add the retained gilts from the production nucleus to the PN pop
    ProdNucPop <-rbind.fill(ProdNucPop,PN_gilts)
    
    print(paste0("PN pop stats: M ",sum(ProdNucPop$sex=="M" & ProdNucPop$fate==1),", F ",sum(ProdNucPop$sex=="F" & ProdNucPop$fate==1)))
    
    # run the ageing for the production nucleus
    ProdNucPop <- ageing(ProdNucPop)
    
    
    #### Multiplier tier ####

    # run cull BEFORE animals are tranferred in and potential parents are selected
    MultPop <- BoarCull(MultPop, AgeDist_M)
    MultPop <- SowCull(MultPop, AgeDist)

    # add in gilts from Prod nucleus, males from SPF B
    MultPop <- rbind.fill(MultPop,Mult_gilts,Mult_boars)
    
    MultPop$merit <- as.numeric(as.character(MultPop$merit))
    MultPop$ID <- as.character(MultPop$ID)
    MultPop$breed <- as.character(MultPop$breed)
    MultPop$age <- as.numeric(as.character(MultPop$age))
    
    
    Mult_males <- subset(MultPop,MultPop$sex=="M" & MultPop$breed=="B" & MultPop$age >= M_ageFirstMate & MultPop$fate==1)
    Mult_dams <- subset(MultPop,MultPop$sex=="F" & MultPop$age >= F_ageFirstMate & MultPop$age %% FarrowInt==rem & MultPop$fate==1)
    
    # create new generation of prod nucleus, top gilts will be retained, and gilts transferred to the multiplier
    newGenMult <- createOffspringPRRS(Mult_males,Mult_dams,LitterSize,paste0("M",g,"_"),indexSD,g)
    pmvec[5] <- mean(as.numeric(as.character(newGenMult$merit)), na.rm = TRUE)
    newGenMult_F <- subset(newGenMult,newGenMult$sex=="F")
    newGenMult_F <- newGenMult_F[order(newGenMult_F$merit,decreasing = TRUE),]

    print(paste0("Mult pop stats: M ",sum(MultPop$sex=="M" & MultPop$fate==1),", F ",sum(MultPop$sex=="F" & MultPop$fate==1)))
    
    if("Mult" %in% genoT){
      # if genotyping in the multiplier tier, select gilts to go to BW based on PRRS resistant first
      BW_gilts <- retainGeno(newGenMult_F,selnF_M_BW)
    } else {
      # otherwise select on merit alone
      BW_gilts <- newGenMult_F[1:selnF_M_BW,]
    }
    
    # run the ageing for the multiplier
    MultPop <- ageing(MultPop)
    
    #### Breeder weaner ####
    # run the cull BEFORE animals are transferred in and potential parents are selected
    BWPop <- BoarCull(BWPop, AgeDist_M)
    BWPop <- SowCull(BWPop, AgeDist)
    
    BWPop <- rbind.fill(BWPop,BW_gilts,BW_boars)
    
    BWPop$age <- as.numeric(as.character(BWPop$age))
    BWPop$merit <- as.numeric(as.character(BWPop$merit))
    BWPop$breed <- as.character(BWPop$breed)
    
    BW_males <- subset(BWPop,BWPop$sex=="M" & BWPop$breed=="T" & BWPop$age >= M_ageFirstMate & BWPop$fate==1)
    BW_dams <- subset(BWPop,BWPop$sex=="F" &  BWPop$breed=="BA" & BWPop$age >= F_ageFirstMate & BWPop$age %% FarrowInt==rem & BWPop$fate==1)
    
    # run the ageing for the breeder weaner
    BWPop <- ageing(BWPop)
    
    newGenBW <- createOffspringPRRS(BW_males,BW_dams,LitterSize,paste0("BW",g,"_"),indexSD,g)
    pmvec[6] <- mean(as.numeric(as.character(newGenBW$merit)), na.rm = TRUE)
    pmvec[7] <- nrow(newGenBW)
    
    print(paste0("BW pop stats: M ",sum(BWPop$sex=="M" & BWPop$fate==1),", F ",sum(BWPop$sex=="F" & BWPop$fate==1)))
    
    # The piglets born in the breeder weaner are commercial slaughter animals - we assume that these would now go into Jaap's code as the pool of pigs 
    CommPop <- rbind(CommPop,newGenBW)
    
    ## Print stats on the number of homozgote resistant pigs in each tier per batch 
    
    print(paste0("Number of homozygote resistant pigs in SPF: ",sum(na.omit(SPFpop$geno=="pp"))))
    print(paste0("Number of homozygote resistant pigs in PN: ",sum(na.omit(ProdNucPop$geno=="pp"))))
    print(paste0("Number of homozygote resistant pigs in Mult: ",sum(na.omit(MultPop$geno=="pp"))))
    nResistant <- sum(na.omit(CommPop$geno=="pp"))
    print(paste0("Number of homozygote resistant pigs in commercial: ",nResistant))
    
    # record the proportion of the population that is carrying the homozygote resistant genotype
    ResistVec <- c(sum(na.omit(SPFpop$geno=="pp")),sum(na.omit(ProdNucPop$geno=="pp")),sum(na.omit(MultPop$geno=="pp")),nResistant,sum(na.omit(newGenBW$geno=="pp")),nrow(CommPop))
    propresvec[1] <- sum(na.omit(SPFpop$geno=="pp"))/nrow(SPFpop)
    propresvec[2] <- sum(na.omit(ProdNucPop$geno=="pp"))/nrow(ProdNucPop)
    propresvec[3] <- sum(na.omit(MultPop$geno=="pp"))/nrow(MultPop)
    propresvec[4] <- sum(na.omit(newGenBW$geno=="pp"))/nrow(newGenBW)
    
    # add the results vectors to the results data frames
    pigletmerit <- rbind(pigletmerit, pmvec)
    edit_numbers<- rbind(edit_numbers, editvec)
    edit_success <- rbind(edit_success,successvec)
    resistant_vec <- rbind(resistant_vec,ResistVec)
    propres_vec <- rbind(propres_vec,propresvec)
  }
  
  
  # Rename the columns for the results data frames
  colnames(pigletmerit) <- c("SPF_A", "SPF_B", "SPF_T", "PN", "Mult", "BW","nCommBorn")
  colnames(edit_numbers) <- c("SPF_A", "SPF_B", "SPF_T", "PN", "Mult", "BW")
  colnames(edit_success) <- c("SPF_A", "SPF_B", "SPF_T", "PN", "Mult", "BW")
  colnames(resistant_vec) <- c("SPF", "PN", "Mult", "BW","BWborn","SizeComm")
  colnames(propres_vec) <- c("SPF%", "PN%", "Mult%", "BWBrn%","spare1","spare2")
  
  resistant_out <- cbind(resistant_vec,propres_vec)
  
  # For a single replicate, output the data to csv 
  outputdat <- cbind(rep(1:timeSteps),pigletmerit,edit_numbers,edit_success,resistant_out) 
  #write.csv(outputdat,paste0("OutputVectors_ep_",SPFeditProp,"_",Sys.Date(),".csv"),row.names = FALSE,quote = FALSE)
  
  SPFpop$popname <- "SPFNUC"
  ProdNucPop$popname <- "PRODNUC"
  MultPop$popname <- "MULT"
  BWPop$popname <- "BW"
  AllPop2 <- rbind.fill(SPFpop,ProdNucPop,MultPop,BWPop)
  
  # Everything <- rbind.fill(AllPop,AllPop2)
  
  # write the burn in pop to file for faster testing - not usually neccessary unless more detailed inspection of the simulated population is desired
  #write.csv(AllPop2,paste0("TestPopulation_ForwardGen",Sys.Date(),".csv"),row.names = FALSE,quote = FALSE)
  #write.csv(Everything,paste0("TestPopulation_ForwardGenandburnin",Sys.Date(),".csv"),row.names = FALSE,quote = FALSE)
  return(list(pmerit=pigletmerit,eNumbers=edit_numbers,eSuccess=edit_success,nRes=resistant_out))
}