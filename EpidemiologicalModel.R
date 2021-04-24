## Epidemiological simulation model for the paper "Gene editing in Farm Animals: A Step Change for Eliminating Epidemics on our Doorstep?" by Petersen et al. 
## This code refers to the "epidemiological simulation model" described in the above paper, and can be used to generate the data for the Figures 1-3 in the main paper, and Figs. S1 and S2
## See example for simulation scenarios at the bottom of the file


library(plotrix)

## Global constants
cTotalPopulation <- 12e6
cMinHerdSize     <- 10 # override functions that try to create herds smaller than this and put to this constant value
cCheckInput      <- TRUE
cScenarioChoices <- c("UniformEdits","RandomRandomEdits","RandomUniformEdits","ConcentrateEditsRandom","ConcentrateEditsRankedDecrease","ConcentrateEditsRankedIncrease")
# this refers to the distribution of genetically resistant animals between the herds. In the paper only results for the following scenarios are presented:
# RandomUniformEdits = "comprehensive"; ConcentrateEditsRandom = either Optimum or Concentrated distribution; depending on whether the nr of animals is flexible or fixed; RandomRandomEdits = "unregulated" 
cDistChoices     <- c("uniform","normal") # refers to distribution of herd size; in the paper results for 'normal' distribution are presented
cVacScenarios    <- c("EditOnly","VacOnly","EditOrVac","EditAndVac")

## Functions


###########################################
# Create datastructure to hold properties of all replicates (private)
createThresholds <- function(nreps){
  Thresholds    <- matrix(nrow=nreps,ncol=4)
  colnames(Thresholds) <- c("rep","MinNoEdits","MinPropEdits", "PropEditedHerds")
  Thresholds[,] = rep(-99, nreps)
  return(Thresholds)
}
#############################################

# Create datastructure to hold properties of all herds (private)
createHerds <- function(HerdCount){
  Herds    <- matrix(nrow=HerdCount,ncol=5)
  colnames(Herds) <- c("Count","EditCount","R0","R0prime","Infected")
  return(Herds)
}

# create herds with uniform herd size 
getUniformHerds <- function(PopSize,HerdCount){
  if (cCheckInput) {
    if (PopSize < cMinHerdSize) {stop("PopSize cannot be smaller than min herd size")}
    if (HerdCount < 1)          {stop("HerdCount cannot be smaller than 1")}
  }
  HerdSize  <- trunc(PopSize/HerdCount,0)
  if (HerdSize < cMinHerdSize){stop("Requested herd size is smaller than min herd size")}
  Herds     <- createHerds(HerdCount = HerdCount)
  Herds[,1] <- rep(HerdSize,HerdCount)
  Herds     <- as.data.frame(Herds)
  return(Herds)
}

# create herds with normally distributed herd size
getNormalDistHerds <- function(PopSize,HerdCount,HerdSizeSD){
  if (cCheckInput) {
    if (PopSize < cMinHerdSize) {stop("PopSize cannot be smaller than min herd size")}
    if (HerdCount < 1)          {stop("HerdCount cannot be smaller than 1")}
    if (HerdSizeSD < 0)         {stop("HerdSizeSD cannot be negative")}
  }  
  # create initial herds
  HerdSizeMean  <- trunc(PopSize/HerdCount,0)
  Herds         <- createHerds(HerdCount = HerdCount)
  Herds[,1]     <- round(rnorm(n=HerdCount,mean=HerdSizeMean,sd=HerdSizeSD))
  Herds[,1][Herds[,1]<cMinHerdSize] <- cMinHerdSize
  # probability of having extactly the total amount of animals is small: correct this
  CreatedPopSize <- sum(Herds[,1])
  if (CreatedPopSize > PopSize) {
    TooMany <- CreatedPopSize-PopSize
    if (sum(Herds[,1][Herds[,1]>cMinHerdSize])<TooMany){stop("Requested population distribution not possible")}
    while (TooMany>0){
      Candidate <- sample(x=(1:nrow(Herds)),size=1)
      if (Herds[Candidate,1] > cMinHerdSize) {
        Herds[Candidate,1] <- Herds[Candidate,1] - 1
        TooMany <- TooMany -1
      }
    }
  }
  
  if (CreatedPopSize < PopSize){
    TooSmallHerds <- sample(x=(1:HerdCount),size=(PopSize-CreatedPopSize),replace = TRUE)
    for (i in 1:length(TooSmallHerds)){Herds[TooSmallHerds[i],1] <- Herds[TooSmallHerds[i],1] + 1}
    # (loop is required as there may be more individuals to add than there are herds, therefor replace = TRUE when sampling)
  }
  
  Herds         <- as.data.frame(Herds)
  return(Herds)
}

# distribute edited animals uniform over all herds to fraction P (from an indefinite supply)
distributeUniformEdits <- function(Herds,P){
  if (cCheckInput) {if ((P < 0) | (P > 1))     {stop("P should be between 0 and 1")}}
  Herds[,"EditCount"] <- Herds[,"Count"]*P
  return(Herds)
}



# Comprehensive distribution scenario: Distribute edits from a limited supply randomly over all herds 
# (function works insensitive to herd sizes)
distributeRandomUniformEdits <- function(Herds,TotalEditCount){
  if (cCheckInput) {if (TotalEditCount < 0)       {stop("TotalEditCount cannot be smaller than 0")}}
  Herds[,"EditCount"] <- rep(0,nrow(Herds))
  TotalEditCount <- min(TotalEditCount,sum(Herds[,"Count"] ))
  PositionCount <- sum(Herds[,"Count"])
  Positions     <- sample(x=1:PositionCount,size=TotalEditCount, replace = FALSE) #create vector length Total Edit Counts of random positions
  Edits         <- c(rep(0,PositionCount))
  Edits[Positions] <- 1 #indicate which animals get edited by assigning a 1 
  Pointer <-1
  for (i in 1:nrow(Herds)){
    Herds[i,"EditCount"] <- sum(Edits[Pointer:(Pointer+Herds[i,"Count"]-1)])
    Pointer <- Pointer + Herds[i,"Count"]
  }
  return(Herds)
}

# Unregulated distribution scenario: Distribute edits from unlimited supply concentrated over a random fraction of the herds, in random counts.  
# Both "randoms" may be overrruled by filling in a required fraction.
distributeRandomRandomEdits <- function(Herds,HerdFraction = "random",P="random"){
  if (cCheckInput) {
    if (!is.numeric(HerdFraction) & !(HerdFraction=="random"))                   {stop("HerdFraction should either be fraction or random")}
    if (is.numeric(HerdFraction) & ((HerdFraction>1) | (P<0)))                   {stop("HerdFraction should be between 0 and 1")}
    if (!is.numeric(P) & !(P=="random"))                                         {stop("P should either be fraction or random")}
    if (is.numeric(P) & ((P>1) | (P<0)))                                         {stop("P should be between 0 and 1")}
  }
  if (HerdFraction=="random") {HerdFraction <- runif(n=1)}
  SelectedHerds <- c(1:nrow(Herds))
  SelectedHerds <- sample(SelectedHerds,size=round(HerdFraction*nrow(Herds)),replace = FALSE)
  Herds[,"EditCount"] <- 0
  if (P=="random") {
    Herds[SelectedHerds,"EditCount"] <- round(Herds[SelectedHerds,"Count"]*runif(n=length(SelectedHerds)))  
  }else{
    Herds[SelectedHerds,"EditCount"] <- round(Herds[SelectedHerds,"Count"]*P)  
  }
   
  return(Herds)
}


# Concentrated and Optimum Distribution scenarios: Distribute edits from limited supply concentrated over a fraction of the herds (random or sorted by size)
# P can either be provided as a "fixed" fraction (for Concentrated Scenario) or as "flex" for (optimum Scenario), meaning that P is adapted to match the R0 of the herd such that R0' <=1. "flex" refers to the Optimum distribution scenario
  distributeConcentrateEdits <- function(Herds,P,TotalEditCount,HerdOrdering,Decreasing){ 
  if (cCheckInput) {
    if (!is.numeric(P) & !(P=="flex"))                     {stop("P should either be fraction or flex")}
    if (is.numeric(P) & ((P>1) | (P<0)))                   {stop("P should be between 0 and 1")}
    if (TotalEditCount < 0)                                {stop("TotalEditCount cannot be smaller than 0")}
    if (!HerdOrdering=="random" & !HerdOrdering=="ranked") {stop("HerdOrdering should either be random or ranked")} #ranked implies that herds receive edits based on their size
  }
  Herds[,"EditCount"] <- rep(0,nrow(Herds))
  TotalEditCount <- min(TotalEditCount,sum(Herds[,"Count"] ))
  if (HerdOrdering=="random") {Ordering <- sample(x=c(1:nrow(Herds)), size=nrow(Herds))}
  if (HerdOrdering=="ranked") {Ordering <- order(x=Herds[,"Count"],decreasing = Decreasing)}
  i=1
  while ((TotalEditCount>0) && (i<=nrow(Herds))){
    if (P=="flex"){Pherd <- max(0, (1-(1/Herds[Ordering[i],"R0"])))}
    else {Pherd <- P}
    #EditCount <- min((trunc(Herds[Ordering[i],"Count"] * Pherd)+1),Herds[Ordering[i],"Count"]) # removed "+1" in trunc(Herds*Pherd) +1 as it added edits to herds where there were none
    EditCount <- min((round(Herds[Ordering[i],"Count"] * Pherd)),Herds[Ordering[i],"Count"])
    if(EditCount ==1) {EditCount = 0} #correct for adding 1 in the truncation above
    if (EditCount<=TotalEditCount) { #i.e edits left to put into that herd
      Herds[Ordering[i],"EditCount"] <- EditCount
      TotalEditCount       <- TotalEditCount - EditCount
    }
    i <- i + 1
  } 
  return(Herds)
}

# Distribute uniform R0 settings over the herds
distributeUniformR0s <- function(Herds,R0){
  if (cCheckInput) {if (R0 <= 0)       {stop("R0 should be higher than 0")}}
  Herds[,"R0"] <- rep(R0,nrow(Herds))
  return(Herds)
}

# Distribute normal distributed R0 settings over the herds
distributeNormalDistR0s <- function(Herds,R0mean,R0sd,R0min,R0max=NULL){
  if (cCheckInput) {
    if (R0mean <= 0)     {stop("R0mean should be higher than 0")}
    if (R0sd   <= 0)     {stop("R0sd should be higher than 0")}
    if (R0min  <= 0)     {stop("R0min should be higher than 0")}
    if (R0max  <= 0)     {stop("R0max should be higher than 0")}
    if (R0min  >= R0max) {stop("R0min should be lower than R0max")}
    if (R0min  >  R0mean){stop("R0min should be <= R0mean")}
    if (R0max  <  R0mean){stop("R0max should be >= R0mean")}
  }
  Herds[,"R0"] <- rnorm(nrow(Herds),mean=R0mean,sd=R0sd)
  Herds[,"R0"][Herds[,"R0"]<R0min] <- R0min
  if (!is.null(R0max)){Herds[,"R0"][Herds[,"R0"]>R0max] <- R0max}
  return(Herds)
}

# calculate R0' from R0 with edits and vaccines, based on eqn [1] in paper 
calcDirectR0prime <- function(R0,P_e,Epsilon_e,P_v,Epsilon_v){
  if (cCheckInput) {
    if (any(R0  <= 0))                          {stop("R0 should be higher than 0")}
    if (any(!is.numeric(P_e)))                        {stop("P_e should be numeric ")}
    if (any(P_e < 0))                           {stop("P_e should be larger or equal to 0")}
    if (any(P_e > 1))                           {stop("P_e should be smaller or equal to 1")}
    if (any(Epsilon_e < 0) | any(Epsilon_e > 1))   {stop("Epsilon_e should be between 0 and 1")}
    if (any(P_v < 0) | any(P_v > 1))               {stop("P_v should be between 0 and 1")}
    if (any(Epsilon_v < 0) | any(Epsilon_v > 1))   {stop("Epsilon_v should be between 0 and 1")}
  }
  R0prime <- (1-(P_e*Epsilon_e)-(1-P_e)*(P_v*Epsilon_v))*R0
  #this assumes that P_e and P_v are correctly specified in CalculateHerdR0primesVac according to the chosen scenarios
  return(R0prime)
}



##new function that includes different vaccine scenarios
calculateHerdR0primesVac <- function(Herds, Epsilon_e,Epsilon_v,P_exp,VacScenario){
  P_e        <- (Herds[,"EditCount"]/Herds[,"Count"])
  P_e[P_e>1] <- 1
  P_v=rep(99,nrow(Herds))
  #vaccination scenarios added by Andrea: specifies the vector P_v for CalcDirectR0prime
  if (VacScenario =="EditOnly")   {  P_v=rep(0,nrow(Herds))}
  if (VacScenario =="VacOnly")   {
    P_e=rep(0,nrow(Herds))
    P_v = rep(1,nrow(Herds))}
  if (VacScenario =="EditOrVac")   {
    for (i in 1:nrow(Herds)) {
      if (Herds[i,"EditCount"] >0){ P_v[i] = 0} #editing is going on in this herd, so P_e is as above and P_v = 0
        else if (Herds[i,"EditCount"] == 0) {P_v[i] = 1} #no editing in herd i, so mass vaccination
    }
  }
  if (VacScenario =="EditAndVac")   {P_v = rep(1,nrow(Herds))} #P_e as specified above
  
  #calculate nr of herds with edits and vaccination
  nEdits <- Herds[,"EditCount"]
  nHerdsEdit <- length(nEdits[nEdits > 0])
  nHerdsEditFract <- nHerdsEdit/nrow(Herds)
  nHerdsEdit2 <- length(P_e[P_e>0])
  if(nHerdsEdit2 != nHerdsEdit) {stop("nHerdsEdit2 != nHerdsEdit")}
  nHerdsVac <- length(P_v[P_v>0])
  nHerdsVacFrac <- nHerdsVac/nrow(Herds)
  
  print(paste("total fraction  of edited herds:",nHerdsEditFract, " total fraction of vaccinated herds: ", nHerdsVacFrac) )
  Herds[,"R0prime"] <- calcDirectR0prime(R0=Herds[,"R0"],P_e= P_e,
                                         Epsilon_e = rep(Epsilon_e,nrow(Herds)),P_v=P_v,Epsilon_v=rep(Epsilon_v,nrow(Herds)))
  # set R0prime of P_exp proportion of randomly chosen herds to zero
  X_uni = runif(nrow(Herds),0,1)  
  for (i in 1:nrow(Herds)){if (X_uni[i]>P_exp) { Herds[i,"R0prime"] <- 0}}
  return(Herds)
}

testPlotDistributions <- function(HerdCount,HerdDist,TotalEditCount,Scenario,VacScenario, R0,R0Dist,P_e,P_exp, PlotSeparate) {
  if (cCheckInput) {
    if (HerdCount<2)                             {stop("There should be min 2 herds")}
    if (HerdCount>cTotalPopulation)              {stop("HerdCount cannot be larger than cTotalPopulation")}
    if (!(HerdDist %in% cDistChoices))           {stop(paste(HerdDist," is not a valid value for HerdDist within this unit. Choose from:",paste(cDistChoices,collapse=", ")))}
    if (TotalEditCount<0)                        {stop("TotalEditCount cannot be smaller than 0")}
    if (!Scenario %in% cScenarioChoices)         {stop(paste(Scenario," is not a valid value for Scenario within this unit. Choose from:",paste(cScenarioChoices,collapse=", ")))}
    if (!VacScenario %in% cVacScenarios)         {stop(paste(VacScenario," is not a valid value for Scenario within this unit. Choose from:",paste(cVacScenarios,collapse=", ")))}
    if (R0<1)                                    {stop("R0 cannot be smaller than 1")}
    if (!(R0Dist %in% cDistChoices))             {stop(paste(R0Dist," is not a valid value for R0Dist within this unit. Choose from:",cDistChoices))}
    if (!is.numeric(P_e) & !(P_e=="flex"))       {stop("P_e should either be fraction or flex")}
    if (is.numeric(P_e) & ((P_e>1) | (P_e<0)))   {stop("P_e should be between 0 and 1")}
  }
  if (HerdDist=="uniform"){Herds <- getUniformHerds(PopSize = cTotalPopulation,HerdCount = HerdCount)}
  if (HerdDist=="normal") {Herds <- getNormalDistHerds(PopSize = cTotalPopulation,HerdCount = HerdCount,HerdSizeSD = 1000)}
  if (R0Dist=="uniform")  {Herds <- distributeUniformR0s(Herds = Herds, R0 = R0)}
  if (R0Dist=="normal")   {Herds <- distributeNormalDistR0s(Herds,R0mean=R0,R0sd=1,R0min=0.1,R0max=13)}
  if (Scenario == "UniformEdits")                  {Herds <- distributeUniformEdits(Herds = Herds, P=TotalEditCount/sum(Herds[,1]))}
  if (Scenario == "RandomRandomEdits")             {Herds <- distributeRandomRandomEdits(Herds = Herds)}
  if (Scenario == "RandomUniformEdits")            {Herds <- distributeRandomUniformEdits(Herds = Herds, TotalEditCount = TotalEditCount)}
  if (Scenario == "ConcentrateEditsRandom")        {Herds <- distributeConcentrateEdits(Herds = Herds, P=P_e,TotalEditCount = TotalEditCount,HerdOrdering="random")}
  if (Scenario == "ConcentrateEditsRankedDecrease"){Herds <- distributeConcentrateEdits(Herds = Herds, P=P_e,TotalEditCount = TotalEditCount,HerdOrdering="ranked",Decreasing = TRUE)}
  if (Scenario == "ConcentrateEditsRankedIncrease"){Herds <- distributeConcentrateEdits(Herds = Herds, P=P_e,TotalEditCount = TotalEditCount,HerdOrdering="ranked",Decreasing = FALSE)}
  
  Herds <- calculateHerdR0primesVac(Herds=Herds, Epsilon_e = Epsilon_e, Epsilon_v = Epsilon_v, P_exp = P_exp, VacScenario = VacScenario)
  if(!PlotSeparate){par(mfrow=c(2,2))}
  hist(x=Herds[,"Count"],                     xlab="Herd sizes"  ,cex.main=0.5,main=paste("#Herds=",HerdCount,"/HerdDist=",HerdDist,"/#Edits"=TotalEditCount,"/Scenario=",Scenario,"/R0Dist=",R0Dist,"/Pe=",P_e,sep="",collapse=""))
  hist(x=Herds[,"EditCount"],                 xlab="# Edits"     ,cex.main=0.5,main=paste("#Herds=",HerdCount,"/HerdDist=",HerdDist,"/#Edits"=TotalEditCount,"/Scenario=",Scenario,"/R0Dist=",R0Dist,"/Pe=",P_e,sep="",collapse=""))
  hist(x=Herds[,"EditCount"]/Herds[,"Count"], xlab="Realized P_e",cex.main=0.5,main=paste("#Herds=",HerdCount,"/HerdDist=",HerdDist,"/#Edits"=TotalEditCount,"/Scenario=",Scenario,"/R0Dist=",R0Dist,"/Pe=",P_e,sep="",collapse=""))
  hist(x=Herds[,"R0"],                        xlab="R0"          ,cex.main=0.5,main=paste("#Herds=",HerdCount,"/HerdDist=",HerdDist,"/#Edits"=TotalEditCount,"/Scenario=",Scenario,"/R0Dist=",R0Dist,"/Pe=",P_e,sep="",collapse=""))
  hist(x=Herds[,"R0prime"],                   xlab="R0prime"     ,cex.main=0.5,main=paste("#Herds=",HerdCount,"/HerdDist=",HerdDist,"/#Edits"=TotalEditCount,"/Scenario=",Scenario,"/R0Dist=",R0Dist,"/Pe=",P_e,sep="",collapse=""))
  #hist(x=Herds[,"Infected"],                  xlab="Infected"    ,cex.main=0.5,main=paste("#Herds=",HerdCount,"/HerdDist=",HerdDist,"/#Edits"=TotalEditCount,"/Scenario=",Scenario,"/R0Dist=",R0Dist,"/Pe=",P_e,sep="",collapse=""))
}

runSimulationSeries <- function(HerdCount,HerdDist,EditCount_choices,Scenario,VacScenario, R0_choices,R0Dist,P_e,Epsilon_e,P_v,Epsilon_v, P_exp) {
  if (cCheckInput) {
    if (HerdCount<2)                                 {stop("There should be min 2 herds")}
    if (HerdCount>cTotalPopulation)                  {stop("HerdCount cannot be larger than cTotalPopulation")}
    if (!(HerdDist %in% cDistChoices))               {stop(paste(HerdDist," is not a valid value for HerdDist within this unit. Choose from:",paste(cDistChoices,collapse=", ")))}
    if (length(EditCount_choices)<1)                 {stop("No EditCount choices")}
    if (any(EditCount_choices<0))                    {stop("EditCount cannot be smaller than 0")}
    if (!Scenario %in% cScenarioChoices)             {stop(paste(Scenario," is not a valid value for Scenario within this unit. Choose from:",paste(cScenarioChoices,collapse=", ")))}
    if (!VacScenario %in% cVacScenarios)             {stop(paste(VacScenario," is not a valid value for Scenario within this unit. Choose from:",paste(cVacScenarios,collapse=", ")))}
    if (length(R0_choices)<1)                        {stop("No R0 choices")}
    if (any(R0_choices<=0))                          {stop("R0 cannot be equal or smaller than 0")}
    if (!(R0Dist %in% cDistChoices))                 {stop(paste(R0Dist," is not a valid value for R0Dist within this unit. Choose from:",paste(cDistChoices,collapse=", ")))}
    if (!is.numeric(P_e) & !(P_e=="flex"))           {stop("P_e should either be fraction or flex")}
    if ((P_e=="flex") && (Scenario %in% c("UniformEdits","RandomUniformEdits")))
                                                     {stop("P_e flex cannot be combined with Scenario UniformEdits or RandomUniformEdits") }
    if (is.numeric(P_e) & ((P_e>1) | (P_e<0)))       {stop("P_e should be between 0 and 1")}
    if ((Epsilon_e>1) | (Epsilon_e<0))               {stop("Epsilon_e should be between 0 and 1")}
    if ((P_v>1) | (P_v<0))                           {stop("P_v should be between 0 and 1")}
    if ((Epsilon_v>1) | (Epsilon_v<0))               {stop("Epsilon_v should be between 0 and 1")}
    if (is.numeric(P_exp) & ((P_exp>1) | (P_exp<0))) {stop("P_exp should be between 0 and 1")}
  } 
  
  if (HerdDist=="uniform"){Herds <- getUniformHerds(PopSize = cTotalPopulation,HerdCount = HerdCount)}
  if (HerdDist=="normal") {Herds <- getNormalDistHerds(PopSize = cTotalPopulation,HerdCount = HerdCount,HerdSizeSD = 1000)}
  
  Results <- matrix(nrow=length(R0_choices),ncol=1+length(EditCount_choices))
  colnames(Results) <- c("R0",EditCount_choices)
  
  for (R0_index in 1:length(R0_choices)){
    R0 <- R0_choices[R0_index]
    if (R0Dist=="uniform")  {Herds <- distributeUniformR0s(Herds = Herds, R0 = R0)}
    if (R0Dist=="normal")   {Herds <- distributeNormalDistR0s(Herds,R0mean=R0,R0sd=1,R0min=0.1,R0max=13)}
    Results[R0_index,1] <- R0
    for (EC_index in 1:length(EditCount_choices)){
      if (Scenario == "UniformEdits")                  {Herds <- distributeUniformEdits(Herds = Herds, P=P_e)}
      if (Scenario == "RandomRandomEdits")             {Herds <- distributeRandomRandomEdits(Herds = Herds)}
      if (Scenario == "RandomUniformEdits")            {Herds <- distributeRandomUniformEdits(Herds = Herds, TotalEditCount = EditCount_choices[EC_index])}
      if (Scenario == "ConcentrateEditsRandom")        {Herds <- distributeConcentrateEdits(Herds = Herds, P=P_e,TotalEditCount = EditCount_choices[EC_index],HerdOrdering="random")}
      if (Scenario == "ConcentrateEditsRankedDecrease"){Herds <- distributeConcentrateEdits(Herds = Herds, P=P_e,TotalEditCount = EditCount_choices[EC_index],HerdOrdering="ranked",Decreasing = TRUE)}
      if (Scenario == "ConcentrateEditsRankedIncrease"){Herds <- distributeConcentrateEdits(Herds = Herds, P=P_e,TotalEditCount = EditCount_choices[EC_index],HerdOrdering="ranked",Decreasing = FALSE)}
      Herds <- calculateHerdR0primesVac(Herds=Herds,Epsilon_e = Epsilon_e,Epsilon_v = Epsilon_v, P_exp = P_exp, VacScenario = VacScenario)
      R0primes <- Herds[,"R0prime"]
      NoOutbreakCount <- length(R0primes[R0primes <= 1])
      NoOutbreakFract <- NoOutbreakCount/nrow(Herds)
      Results[R0_index,EC_index+1] <- NoOutbreakFract
    }
    print(paste("Completed", R0_index,"/",length(R0_choices)))
  }
  return(Results)
} 



#function to generate the data for figures 2-3, and S1 and S2 in the paper
# this function uses Newton-Raphson solver to find required number of edits to eradicate the disease, i.e. for a stipulated target fraction of herds to have R0'<1
findMinEdits4HerdImmunityFraction <- function(HerdCount,HerdDist,Scenario,VacScenario, R0,R0Dist,
                                              P_e,Epsilon_e,Epsilon_v,P_exp, TargetFraction,
                                              Prior,Tolerance,ConvThres, Verbose) {
  if (cCheckInput) {
    if (HerdCount<2)                                {stop("There should be min 2 herds")}
    if (HerdCount>cTotalPopulation)                 {stop("HerdCount cannot be larger than cTotalPopulation")}
    if (!(HerdDist %in% cDistChoices))              {stop(paste(HerdDist," is not a valid value for HerdDist within this unit. Choose from:",paste(cDistChoices,collapse=", ")))}
    if (!Scenario %in% cScenarioChoices)            {stop(paste(Scenario," is not a valid value for Scenario within this unit. Choose from:",paste(cScenarioChoices,collapse=", ")))}
    if (!VacScenario %in% cVacScenarios)            {stop(paste(VacScenario," is not a valid value for Scenario within this unit. Choose from:",paste(cVacScenarios,collapse=", ")))}
    if (R0<=0)                                      {stop("R0 cannot be equal or smaller than 0")}
    if (!(R0Dist %in% cDistChoices))                {stop(paste(R0Dist," is not a valid value for R0Dist within this unit. Choose from:",paste(cDistChoices,collapse=", ")))}
    if (!is.numeric(P_e) & !(P_e=="flex"))          {stop("P_e should either be fraction or flex")}
    if ((P_e=="flex") && (Scenario %in% c("UniformEdits"))) {stop("P_e flex cannot be combined with Scenario UniformEdits") }
    if (is.numeric(P_e) & ((P_e>1) | (P_e<0)))      {stop("P_e should be between 0 and 1")}
    if ((Epsilon_e>1) | (Epsilon_e<0))              {stop("Epsilon_e should be between 0 and 1")}
    if ((Epsilon_v>1) | (Epsilon_v<0))              {stop("Epsilon_v should be between 0 and 1")}
    if ((TargetFraction<0) | (TargetFraction>1))    {stop("TargetFraction should be between 0 and 1")}
    if ((Prior<0) | (Prior>1))                      {stop("Prior should be between 0 and 1")}
    if (Tolerance<=0)                               {stop("Tolerance should be > 0")}
    if (ConvThres<=10)                               {stop("Convergence Threshold should be > 10")}
    if (!is.logical(Verbose))                       {stop("Verbose should be logical")}
  } 
  if (Scenario == "UniformEdits") {stop("Function not applicable to Scenario UniformEdits")}
  if(VacScenario == "VacOnly") {stop("Function not applicable to scenario VacOnly = P_e = 0")}
  
  if (HerdDist=="uniform"){Herds <- getUniformHerds(PopSize = cTotalPopulation,HerdCount = HerdCount)}
  if (HerdDist=="normal") {Herds <- getNormalDistHerds(PopSize = cTotalPopulation,HerdCount = HerdCount,HerdSizeSD = 1000)}
  
  if (R0Dist=="uniform")  {Herds <- distributeUniformR0s(Herds = Herds, R0 = R0)}
  if (R0Dist=="normal")   {Herds <- distributeNormalDistR0s(Herds,R0mean=R0,R0sd=1,R0min=0.1,R0max=13)}
  
  
  #set initial values and constants
  TimeOut         <- 30
  Round           <- 0
  Convergence     <- FALSE
  TargetReached   <- FALSE
  seed            <- Prior
  Tolerance       <- Tolerance
  ConvThres       <- ConvThres
  totalmin        <- 0
  totalmax        <- cTotalPopulation
  newmin          <- totalmin
  newmax          <- totalmax
  TempSol         <- cTotalPopulation
  minEdits        <- cTotalPopulation
  TotalEditCount  <- newmin+seed*(newmax - newmin)
  NoOutbreakFract <- 0 #Added by Andrea to stop error message
  delta           <- newmax - newmin

  results    <- matrix(nrow=1,ncol=3)
  colnames(results) <- c("Convergence","MinNoEdits","nEditedHerds")
  results[1,1] <- FALSE
  
  ##check convergence - change convergence criterion when using fractions
  if (delta <= ConvThres) {Convergence <- TRUE 
  stop("ConvThreshold set too large")}
  
  while ((!Convergence) && (Round<=TimeOut)){
    ##Distribute the number of edited individuals into the herds
    if (Verbose) {print(paste("Iteration",Round,": Current fraction=",NoOutbreakFract," > Set edits=",TotalEditCount))}
    ##distribute the edits across the herds as specified by editing scenario
    if (Scenario == "RandomUniformEdits")            {Herds <- distributeRandomUniformEdits(Herds = Herds, TotalEditCount = TotalEditCount)}
    if (Scenario == "RandomRandomEdits")             {Herds <- distributeRandomRandomEdits(Herds = Herds)}
    if (Scenario == "ConcentrateEditsRandom")        {Herds <- distributeConcentrateEdits(Herds = Herds, P=P_e,TotalEditCount = TotalEditCount,HerdOrdering="random")}
    if (Scenario == "ConcentrateEditsRankedDecrease"){Herds <- distributeConcentrateEdits(Herds = Herds, P=P_e,TotalEditCount = TotalEditCount,HerdOrdering="ranked",Decreasing = TRUE)}
    if (Scenario == "ConcentrateEditsRankedIncrease"){Herds <- distributeConcentrateEdits(Herds = Herds, P=P_e,TotalEditCount = TotalEditCount,HerdOrdering="ranked",Decreasing = FALSE)}
    
    ##Calculate the proportion of herds with R0prime <=1
    Herds <- calculateHerdR0primesVac(Herds=Herds, Epsilon_e=Epsilon_e,Epsilon_v=Epsilon_v, P_exp = P_exp, VacScenario =VacScenario)
   
    R0primes <- Herds[,"R0prime"]
    NoOutbreakCount <- length(R0primes[R0primes <= 1])
    NoOutbreakFract <- NoOutbreakCount/nrow(Herds)
    
    ##check if TotalEditCount is within [Target-Tol, Target+Tol]
    if ((NoOutbreakFract >= (TargetFraction-Tolerance)) && (NoOutbreakFract <= (TargetFraction+Tolerance))){ # solution found!
      TargetReached <- TRUE
      TempSol = TotalEditCount
      newmax = TotalEditCount
    }
    else if (NoOutbreakFract>(TargetFraction+Tolerance)) {#too many herds experience no outbreak; reduce nr
      newmax = TotalEditCount
      TargetReached = FALSE
    }
    else if (NoOutbreakFract < (TargetFraction - Tolerance)){
      newmin = TotalEditCount
      TargetReached = FALSE}
    
    ##check convergence
    delta = newmax - newmin
    ##check convergence - change convergence criterion when using fractions
    if (delta<= ConvThres) {Convergence <- TRUE}
    
    ## update solution
    #if (!TargetReached){
      TotalEditCount = newmin+0.5*(newmax-newmin)
      #}
    Round <- Round + 1
  } #while
  
  ##Set final values
  if(TargetReached){minEdits = TotalEditCount
    if (Verbose) {print(paste("Convergence: Final fraction of herds with R0'<=1 =",NoOutbreakFract,"for # edits=",minEdits))}
  }
  else {minEdits <- TempSol
    if (Verbose) {print("No convergence")}
  }
  #return(round(minEdits))
  results[1,1] <- TargetReached 
  results[1,2] <- round(minEdits)
  nHerdsEdit <- length(Herds[,"EditCount"][Herds[,"EditCount"] >0])
  results[1,3] <- nHerdsEdit
  return(results)  
    
}


# functions to calculate averages and SDs over multiple simulation rounds
calcResultsAverages <- function(Results,RepetitionCount){
  if (RepetitionCount<1) {stop("Cannot calculate average for fewer than 1 repetition")}
  NetResultCount <- nrow(Results)/RepetitionCount
  if (!(trunc(NetResultCount)==NetResultCount)){stop("NRow of Results is not a multiple of RepetitionCount") }
  NetResultCount <- nrow(Results)/RepetitionCount
  Av_fracts <- matrix(nrow=NetResultCount,ncol=ncol(Results)-1)
  colnames(Av_fracts) <- colnames(Results)[2:ncol(Results)]
  rownames(Av_fracts) <- Results[1:NetResultCount,1]
  for (i in 1:NetResultCount){
    Indices <- seq(from=i,to=nrow(Results),by=NetResultCount)
    for (j in 2:ncol(Results)){
      Av_fracts[i,j-1] <- mean(as.numeric(Results[Indices,j]))
    }
  }
  return(Av_fracts)
}

calcResultsSDs     <- function(Results,RepetitionCount){
  if (RepetitionCount<2) {stop("Cannot calculate SD for fewer than 2 repetitions")}
  NetResultCount <- nrow(Results)/RepetitionCount
  if (!(trunc(NetResultCount)==NetResultCount)){stop("NRow of Results is not a multiple of RepetitionCount") }
  SD_fracts <- matrix(nrow=NetResultCount,ncol=ncol(Results)-1)
  colnames(SD_fracts) <- colnames(Results)[2:ncol(Results)]
  rownames(SD_fracts) <- Results[1:NetResultCount,1]
  for (i in 1:NetResultCount){
    Indices <- seq(from=i,to=nrow(Results),by=NetResultCount)
    for (j in 2:ncol(Results)){
      SD_fracts[i,j-1] <- sd(as.numeric(Results[Indices,j]))
    }
  }
  return(SD_fracts)
}

# plotting summary over multiple simulation rounds  
linePlotAverageResults <- function(ResultsAverages,ResultsSDs=null,xlab="x",ylab="y",
                                   ShowLegend=TRUE, LegendPrefix="", LegendXpos=0.8, 
                                   InverseYs=FALSE, ByCols=TRUE,ylim=c(0,1)){
  if (!ByCols){
    ResultsAverages <- t(ResultsAverages)
    if (!is.null(ResultsSDs)){ResultsSDs <- t(ResultsSDs)}
  }
  
  cl_line <- rainbow(n=ncol(ResultsAverages),alpha=1)
  cl_CI   <- rainbow(n=ncol(ResultsAverages),alpha=0.1)
  
  # Averages
  Xs <- as.double(rownames(ResultsAverages))
  if (InverseYs){
    Ys <- 1- ResultsAverages
  } else {
    Ys <- ResultsAverages
  }  
  matplot(x=Xs,y=Ys,type = "l",xlab=xlab,ylab=ylab,
          lwd = 1, lty = 1, bty = "l", col = cl_line,ylim=ylim)
  
  # SDs (if available)
  if (!is.null(ResultsSDs)){
    for (i in 1:ncol(Ys)){dispersion(x=Xs,y=Ys[,i],ulim=ResultsSDs[,i],llim=ResultsSDs[,i],type="l",fill=cl_CI[i]) }
  }
  
  # Legend
  if (ShowLegend){
    X_legend <- max(as.double(row.names(ResultsAverages)))*LegendXpos 
    Y_legend <- 1
    legend(x=X_legend, y=Y_legend, paste(LegendPrefix,as.double(colnames(ResultsAverages))), pch = 1,cex=0.5, col = cl_line, bty = "n")
  }
  
}



# I/O fuctions for exchange of results among multiple R terminal instances
# (for parallel runs)
saveRawResults <- function(Results,FileName){
  write.csv(Results,FileName)
}

appendRawResults <- function(Results,FileName){
  AppendResults <- read.csv(FileName)
  AppendResults <- AppendResults[,2:ncol(AppendResults)]
  colnames(AppendResults) <- colnames(Results)
  if (is.null(Results)){
    Results <- AppendResults
  }else{
    Results <- rbind(Results,AppendResults)
  }  
  return(Results)
}




######################################################################################################

##EXECUTION EXAMPLE To Produce results for the figures in the paper

#######################################################################################################
# example parameter set:
HerdCount <- 5000
R0        <- 1.5 #mean R0 value of herds without control

Epsilon_e <- 1 #editing efficiency, i.e. proportion of edits that are resistant
Epsilon_v <- 0.7 #vaccine efficiency
TotalEditCount <- 6e6 #nr of edits available (in principle)
seed <- 0.5
P_exp  <- 1.0 #exposure probability


# choose distribution of herd size and R0
#HerdDist <- "uniform"
HerdDist <- "normal"
#R0Dist   <- "uniform"
R0Dist   <- "normal"

# Choose one Scenario representing the distribution of genetically resistant animals between herds
#Scenario <- "UniformEdits" #not represented in paper; all herds receive equal amount of edits
#Scenario <- "RandomUniformEdits"# = Comprehensive distribution scenario in paper
#Scenario <- "RandomRandomEdits" # = Unregulated distribution scenario
Scenario <- "ConcentrateEditsRandom" # = optimum or Concentrated distribution scenario. Choose P_e = 'flex' for optimum below, or a fixed value of P_e otherwise
#Scenario <- "ConcentrateEditsRankedDecrease"
#Scenario <- "ConcentrateEditsRankedIncrease"

#choose P_e = proportion of edits received by each herd that receives edits
#P_e       <- 1.0 - (1.0/(R0+1.645)) #expected critical prop. edits required for R0prime_avg<1; assumes mean R0 is known or can be estimated
P_e       <- 0.5 # fixed proportion of edits in each of the herds that get edits
#P_e 	<- "flex" #choose 'flex' for optimum distribution scenario (i.e. each herd receives the required nr of edits to achieve R0'<1)

#choose one vaccination scenario
#VacScenario <- "EditOnly"
#VacScenario <- "VacOnly"
#VacScenario <- "EditOrVac"
VacScenario <- "EditAndVac"
P_v       <- 1 # proportion vaccinated in each herd

# Test plot distributions
#######################################################################################################
testPlotDistributions(HerdCount=HerdCount,HerdDist=HerdDist,TotalEditCount=TotalEditCount,Scenario=Scenario,VacScenario=VacScenario,
                      R0=R0, R0Dist=R0Dist, P_e=P_e, P_exp = P_exp, PlotSeparate=FALSE) 
#######################################################################################################


##EXECUTION EXAMPLE To Produce results for FOR FIGURE 1
# Simulation rounds, looping over R0 and Pe

RepetitionCount <- 10 # number of repetitions per parameter setting 

# Choose one of the range options for R0
# Regular run
R0_choices <- rep(seq(from=1,to=10,by=1),RepetitionCount)
# R0 high resolution run (for R0<5)   (slow)
R0_choices <- rep(c(seq(from=1,to=4.9,by=0.1),seq(from=5,to=10,by=1)),RepetitionCount)
# R0 high resolution run (full range) (even slower)
R0_choices <- rep(seq(from=1,to=10,by=0.1),RepetitionCount)

# Set the range of amounts of edits to be tested
EditCount_choices  <- seq(from = 1.0e6, to = 1.2e7, by = 1.0e6)

# run simulations
Results <- runSimulationSeries(HerdCount=HerdCount,HerdDist=HerdDist,EditCount_choices = EditCount_choices, Scenario=Scenario, VacScenario = VacScenario,
                               R0_choices=R0_choices,R0Dist=R0Dist, P_e=P_e,Epsilon_e=Epsilon_e,P_v=P_v,Epsilon_v=Epsilon_v, P_exp = P_exp)


# calc averages over the repetitions
Av_fracts <- calcResultsAverages(Results=Results,RepetitionCount = RepetitionCount)
SD_fracts <- calcResultsSDs     (Results=Results,RepetitionCount = RepetitionCount)
# plot average results
linePlotAverageResults(ResultsAverages = Av_fracts,ResultsSDs = SD_fracts,
                       ShowLegend = TRUE,xlab="R0 mean",ylab="Fraction of no outbreak herds")
#######################################################################################################


##EXECUTION EXAMPLE To Produce results for FOR FIGURES 2, 3
# Get estimates for "How many edits are required to eradicate the disease, i.e. get TotalFreaction of the herds immune?"

#Set up replicates and results matrix
nreps <- 50
Thresholds = createThresholds(nreps)

i <- 1
while (i <= nreps) # wrapper to collect the required number of succesfull attempts. May run forever...
{

  print(paste("(re-)commencing repetition = ",i))
  x_out <- 0
  seed = runif(1,0,1)
  
  x_out <- findMinEdits4HerdImmunityFraction(HerdCount=HerdCount,HerdDist=HerdDist,Scenario=Scenario,VacScenario = VacScenario, R0=R0,
                                             R0Dist=R0Dist,P_e=P_e,Epsilon_e=Epsilon_e,Epsilon_v=Epsilon_v,P_exp,
                                             TargetFraction=0.995,Prior=seed,Tolerance=0.005, ConvThres = 1000, Verbose=TRUE)
  
  if (x_out[1,1]) {
    Thresholds[i,"rep"] = i
    Thresholds[i, "MinNoEdits"] = x_out[1,2]
    Thresholds[i, "MinPropEdits"] = x_out[1,2]/cTotalPopulation
    Thresholds[i, "PropEditedHerds"]= x_out[1,3]/HerdCount
    i <- i + 1
  } 
}
AvgThreshold = calcResultsAverages(Thresholds, nreps) 
SDThreshold = calcResultsSDs(Thresholds, nreps)
SEThreshold =SDThreshold / sqrt(nreps)
print(paste("AvgThreshold -nr", "SEThreshold - nr" ,round(AvgThreshold[1]),round(SEThreshold[1])))
print(paste("AvgThreshold for propEdits and SE",AvgThreshold[2],SEThreshold[2], "AvgPropEditedHerds and SE",AvgThreshold[3],SEThreshold[3]))
print(paste("P_e",P_e))
#######################################################################################################

