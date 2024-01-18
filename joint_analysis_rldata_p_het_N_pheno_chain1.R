#Script_joint_model
library(nimble)
library(nimbleEcology)


#Functions from nimble forums to get and set model and MCMC state and variables (to restart MCMC chain where it left off after saving)
getStateVariableNames <- function(samplerDef) {
  resetMethod <- body(samplerDef$reset)
  stateVars <- character()
  if(resetMethod[[1]] != '{') stop('something wrong')
  numLines <- length(resetMethod)
  for(i in 1:numLines) {
    if(i == 1) next
    thisLine <- resetMethod[[i]]
    if(thisLine[[1]] == '<<-') {
      LHS <- thisLine[[2]]
      if(!is.name(LHS)) stop('haven\'t dealt with non-name-LHS case yet')
      stateVars <- c(stateVars, as.character(LHS))
    }
    if('my_calcAdaptationFactor' %in% all.names(thisLine)) {
      stateVars <- c(stateVars, 'my_calcAdaptationFactor')
    }
  }
  setupMethod <- body(samplerDef$setup)
  if('empirSamp' %in% all.names(setupMethod)) stateVars <- c(stateVars, 'empirSamp')
  return(stateVars)
}

getModelState <- function(model) {
  modelVarNames <- model$getVarNames()
  modelVarValuesList <- vector('list', length(modelVarNames))
  names(modelVarValuesList) <- modelVarNames
  for(var in modelVarNames) {
    modelVarValuesList[[var]] <- model[[var]]
  }
  return(modelVarValuesList)
}

getMCMCstate <- function(conf, mcmc) {
  stateVarNamesList <- vector('list', length(conf$samplerConfs))
  mcmcStateValuesList <- vector('list', length(conf$samplerConfs))
  for(i in seq_along(conf$samplerConfs)) {
    samplerDef <- conf$getSamplerDefinition(i)
    theseStateNames <- getStateVariableNames(samplerDef)
    theseStateValuesList <- vector('list', length(theseStateNames))
    names(theseStateValuesList) <- theseStateNames
    for(j in seq_along(theseStateNames)) {
      if(is.nf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$timesAdapted,
                                            gamma1 = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$gamma1)
        } else
          theseStateValuesList[[j]] <- mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]
      }
      if(is.Cnf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted'),
                                            gamma1 = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1'))
        } else
          theseStateValuesList[[j]] <- valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], theseStateNames[j])
      }
    }
    mcmcStateValuesList[[i]] <- theseStateValuesList
  }
  return(mcmcStateValuesList)
}

#MCMC parmeters
n.thin<-3
n.thin2<-30
n.iter = 1111
n.burnin = 0

#Some useful functions
#Function for 'not in'
'%!in%' <- function(x,y)!('%in%'(x,y))

#Function for calculating difference in number of lines between 2 dataframes
diff<-function(x,y){nrow(x)-nrow(y)}

#Number of unique values in a vector
lu<-function(x){length(unique(x))}

#Load dataset
all.data<-readRDS("./capt_hist_AY_1990_2019_measures_pheno.RDS")

#Annual date of peak nitrogen in plants
pheno<-all.data$N.dat

# Compute time of first capture for each individual
get.first <- function(x) min(which(x!=0))
f.all <- apply(all.data$capt.hist.final[2:31], 1, get.first)

#Ientify ID of birds marked on last occasion and remove them as they are not informative for survival analysis
marked.last.occ<-which(f.all==30)
goose.data<-all.data$capt.hist.final[-marked.last.occ,]

#Sample data to test if model runs on small # of birds (must un-comment the following 4 lines, don't forget to comment back before sending to cluster)
#####################################################################################################################################
# goose.wbtgs<-goose.data[which(goose.data$WT==1),][seq(1,2312,length.out=1000),]
# table(goose.wbtgs$year)
# 
# goose.data<-goose.data[seq(1,nrow(goose.data), length.out=1000),]
# 
# goose.data<-rbind(goose.wbtgs,goose.data)
#####################################################################################################################################


#Capture history
CH<-goose.data[,2:31]

#Determine occasion of marking for each individual
get.first <- function(x) min(which(x!=0))
f <- as.vector(apply(CH, 1, get.first))


#Compute number of individuals
nind<-dim(goose.data)[1]

#Define time indices for each period with different hunting regulations
period1<-1:8     #Historical hunting regulations (1990-1998)
period2<-9:18    #Special measures in Canada only (1998-2008)
period3<-19:30   #Special measures in Canada and the USA (2009-2019)

#Determine sex of each individual (1 male; 2 female) from group number
sex<-ifelse(goose.data$group%in%c(1,3),1,2)
table(sex) #should be more males than females because ~2/3rd of adult females marked with collars and excluded

#Determine # of occasions
n.occasions <- dim(CH)[2]


#Get age class of individuals (adult vs. juv)
age<-ifelse(goose.data$ageB=='ad',2,1)

#Get lineID of juveniles and adults
juv <- which(age==1)
ad  <- which(age==2)

#Extract mismatch value for juveniles (those with known and unknown mismatch values)
ind.mismatch<-goose.data$mismatch[which(age==1)]

#Extract position of individuals with known mismatch value (known age in days since banding, i.e. those marked with webtags)
ID.ms.kwn <- which(goose.data$WT==1)

#Extract position of individuals without known mismatch value ( i.e. juveniles not marked with webtags and adults)
ID.ms.unk<-c(which(goose.data$WT==0),ad)

# Recode CH matrix: note, a 0 is not allowed!
# 1 = alive and in study are, 2 = recovered dead, 3 = not seen or recovered
rCH <- CH # Recoded CH
rCH[rCH==0] <- 3

#Calculate annual average mismatch from known age juveniles
average_yearly_mismatch<-rep(NA, times=max(unique(f)))

for(t in 1:max(f)){
 
  ID_yr_t<-which(f[ID.ms.kwn]==t)
  
  average_yearly_mismatch[t]<-mean(ind.mismatch[ID.ms.kwn][ID_yr_t])
   
}

#No webtags in 1990, and no primaries measured, pretty much impossible to estimate juvenile survival in that year. Just dont make predictions for 1990 for juveniles (but need a value for the model to run, use mean) 
average_yearly_mismatch[1]<-mean(average_yearly_mismatch[2:max(unique(f))])


#Discretize mismatch values to speed up computation time
#Define a mean and sd that are generally in the vicinity of real values to do some scaling
scale.ms.mean<-10.9
scale.ms.sd<-5.9

#Create some potential integer mismatch values (much larger range than observed to accomodate extreme values drawn during MCMC)
pot.ms.vals<-seq(from=-35, to=75)

Scale these potential values with mean and sd defined above
scaled.pot.ms<-(pot.ms.vals-scale.ms.mean)/scale.ms.sd

#Create an index from 1 to length(pot.ms.vals) for each value
index.ms.values<-pot.ms.vals-min(pot.ms.vals)+1

#Create vector to classify individual mismatch values (ranging from -35 to 75) into bins of ~3 days (larger bins at start and end, 3 days in middle)
ms.class.vec<-rep(NA, times=length(scaled.pot.ms))

ms.class.vec[1:30]<-1 #all these values are not observed and this bin will only contain a handful of draws from MCMC so reduce to 1 bin
range.possible.mismatch<-c(-5:33) #observed mismatch values

#Classify potential mismatch values in discrete bins
df.range.possible.mismatch<- data.frame(real.values=range.possible.mismatch,
                                        reclass=c(2,2,2,2,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10))
ms.class.vec[which(pot.ms.vals %in% range.possible.mismatch)]<-df.range.possible.mismatch$reclass
ms.class.vec[(which(pot.ms.vals == max(range.possible.mismatch))+1):length(ms.class.vec)]<-11


#Compute a scaled ms value for each mismatch bin using mean value in bin, and equivalent scaling from dataset
info.bin.ms<-data.frame(real.ms=pot.ms.vals,scaled_ms=scaled.pot.ms, reclass.val=ms.class.vec)

average.bin.ms.scaled<-rep(NA, times=length(unique(ms.class.vec)))

for(i in unique(ms.class.vec)){
  average.bin.ms.scaled[i]<-mean(info.bin.ms$scaled_ms[which(info.bin.ms$reclass.val==i)])
}

# mean of mcmc samples (from tests not in this script) for mismatch values >33 = 36, so we use scaled value of mismatch = 36 to plug in for all mismatches >33
# Similarly, we use mismatch for -7 for all mismatch values < -5
average.bin.ms.scaled[c(1,11)]<-c(-3.10,3.77)
# check<-data.frame(real.ms=pot.ms.vals,scaled_ms=scaled.pot.ms, reclass.val=ms.class.vec)
# check.samples<-table(chain1_run1_age[,1:1300])
# check.samples
# mean(rep(65:77, times=check.samples[43:55]))


#Make sure all potential values correspond to an index
length(index.ms.values)==length(scaled.pot.ms)

#get the range of mismatch values over which to predict survival
range.mismatch.yr<-seq(from=(min(ind.mismatch[ID.ms.kwn])-scale.ms.mean)/scale.ms.sd, 
                       to=(max(ind.mismatch[ID.ms.kwn])-scale.ms.mean)/scale.ms.sd, 
                       length.out=100)

mean(pheno$peakN)
sd(pheno$peakN)

#Get range of mismatch for early, average, and late phenology years to see what predictions look like when considering phenology
pheno$peakN_scaled<-((pheno$peakN-mean(pheno$peakN))/sd(pheno$peakN))


range.mismatch.pheno.early<-seq(from=(range(goose.data[which(goose.data$peakN<176),]$mismatch, na.rm=T)[1]-scale.ms.mean)/scale.ms.sd, 
                            to=(range(goose.data[which(goose.data$peakN<176),]$mismatch, na.rm=T)[2]-scale.ms.mean)/scale.ms.sd, 
                            length.out=15)

range.mismatch.pheno.med<-seq(from=(range(goose.data[which(goose.data$peakN%in%c(176:181)),]$mismatch, na.rm=T)[1]-scale.ms.mean)/scale.ms.sd, 
                            to=(range(goose.data[which(goose.data$peakN%in%c(176:181)),]$mismatch, na.rm=T)[2]-scale.ms.mean)/scale.ms.sd, 
                            length.out=15)

range.mismatch.pheno.late<-seq(from=(range(goose.data[which(goose.data$peakN>181),]$mismatch, na.rm=T)[1]-scale.ms.mean)/scale.ms.sd, 
                            to=(range(goose.data[which(goose.data$peakN>181),]$mismatch, na.rm=T)[2]-scale.ms.mean)/scale.ms.sd, 
                            length.out=15)

scaled.pN.early=(mean(goose.data[which(goose.data$peakN<176),]$peakN)-mean(pheno$peakN))/sd(pheno$peakN)
scaled.pN.med=(mean(goose.data[which(goose.data$peakN%in%c(176:181)),]$peakN)-mean(pheno$peakN))/sd(pheno$peakN)
scaled.pN.late=(mean(goose.data[which(goose.data$peakN>181),]$peakN)-mean(pheno$peakN))/sd(pheno$peakN)



#Index for gamma matrix, where all adult individuals share the same transition matrix slice for a given time (no individual covariates for adults, this savec a lot of computation time for building the model graph)
gi.ad<-length(unique(ms.class.vec))+1

#Data, monitors and constants

#Indices for known and unknown ages of juveniles ONLY (for age model, these vectors differ from ID.ms.kwn and ID.ms.unk)
b.age.kwn.vec<-which(goose.data$WT==1)
b.age.unk.vec<-which(goose.data$WT==0)

#Random sample of ages to track during MCMC for checks
age.keep<-sort(sample(b.age.unk.vec, size=round(min(100,length(b.age.unk.vec)/2)), replace=F))

#Data
known_ages <- goose.data$age.at.B[b.age.kwn.vec]


#Length of primary for juveniles
prim9_lghts <- goose.data$prim9.B[juv]
#Scale it
scaled_prim9_lghts <- as.vector(scale(prim9_lghts))

#Extract mismatch value for juveniles of known age
mismatch.known<-ind.mismatch[ID.ms.kwn]


#List of year of marking for within model mean mismath calculation
juv.list<-list()

for(i in 1:(n.occasions-1)){
  
  juv.list[[i]]<-which(f[juv]==i)
  
}

#Generate matrix with ID of juveniles marked each year
juv.f.indices<-array(dim=c(max(table(goose.data[which(goose.data$ageB=='juv'),]$year)),n.occasions-1))
n.juv.marked.yr<-rep(NA, times=n.occasions-1)
for(t in 2:(n.occasions-1)){
  
  juv.f.indices[1:length(juv.list[[t]]),t]<-juv.list[[t]]
  n.juv.marked.yr[t]<-length(juv.list[[t]])
  
}



#Pack data and constants
#Data
my.data.joint <- list(age.at.B = known_ages,
                      prim9 = as.vector(scale(prim9_lghts)),
                      mismatch.known=mismatch.known,
                      peakN=pheno$peakN,
                      scaled.peakN=pheno$peakN_scaled,
                      scaled.pN.early=scaled.pN.early,
                      scaled.pN.med=scaled.pN.med,
                      scaled.pN.late=scaled.pN.late,
                      julian_date_B_unk = goose.data$julian.date.B[b.age.unk.vec],
                      ms.class.vec = ms.class.vec,     #Index to reclassify mismatch values into 11 categories
                      y=rCH)

#Constants
my.constants.joint <- list(
  #### AGE MODEL ##### #
  n.wbtg=length(b.age.kwn.vec),
  n.not.wbtg=length(b.age.unk.vec),
  b.age.kwn.vec=b.age.kwn.vec,
  b.age.unk.vec=b.age.unk.vec,
  
  ##############################################################################################################################
  #Voir avec Roger comment on gère le fait qu'on a pas de webtags et donc pas possible d'estimer un effet de l'année pour 1990...
  n.yrs = (length(pheno$year)-1),
  ##############################################################################################################################
  
  
  age.keep = age.keep,                           #ID of juveniles whose ages we want to monitor.. was used for simulations.. may not be worth keeping anymore....
  yr.unk = f[b.age.unk.vec],
  yr.kwn = f[b.age.kwn.vec],
  
  ##### SURVIVAL MODEL ##### #
  #Numbers to feed to for loops of different lengths
  n.phunt=2,                                     # number of periods with hunting regulations changes (excluding period #1)
  end.occasion.p1=max(period1),                  # occasion précédant l'implantation des mesures de gestion au Québec en 1999
  end.occasion.p2=max(period2),                  # occasion précédant l'implantation des mesures de gestion aux USA en 2009
  n.juv=length(juv),                             #Total number of juveniles
  #n.ad=length(ad),                              #Total number of adults
  n.age.save = length(age.keep),                 #Number of juveniles for which we want to monitor the estimated age
  n.ms.class = length(unique(ms.class.vec)),     #Number of classes for mismatch effect estimation
  
  #Time indices
  f = f,                                  #occasion of marking (first capture) for all individuals
  n.occasions = dim(rCH)[2],              #Number of capture occasions
  hunt=c(rep(0, times=length(period1)),
         rep(1, times=length(period2)),
         rep(2, times=length(period3))),  # index du beta p.hunt en vigueur à l'occasion en cours (vaut 1 partout pour modèle test, et 1 ou 2 pour modèle à 3 périodes), longueur = n.occasions vaut 0 et ensuite 1 (et ensuite 2 pour vrai modèle)
  #N=N,                                    #Number of juveniles or adults marked at each ocasion (used to compute mean annual juvenile survival)
  
  #ID indices
  nind = dim(rCH)[1],                              #Total number of individuals
  gi.ad=gi.ad,                                     #Index for transition matrix (all adults share the matrix slice indexed at n.juv+1)
  #Sj.pred.index=Sj.pred.index,                    #Index to select only juveniles marked in year t when computing mean annual survival of juveniles
  average_yearly_mismatch=average_yearly_mismatch, #comment this line if using the model with annual mismatch only
  juv.f.indices=juv.f.indices,
  n.juv.marked.yr=n.juv.marked.yr,
  
  #Other constants:
  min.ms = min(pot.ms.vals),               #Feed smallest mismatch to model to be sure to have enough classes
  mismatch.scaled = average.bin.ms.scaled, #Feed scaled mismatch values to put into the appropriate transition matirx slices
  scale.ms.mean=scale.ms.mean,             #Feed mismatch scaling parameters to model
  scale.ms.sd=scale.ms.sd,
  age=age, 
  sex=sex,
  n.sex=2,
  range.mismatch=range.mismatch.yr, 
  range.mismatch.pheno.early=range.mismatch.pheno.early,
  range.mismatch.pheno.med=range.mismatch.pheno.med,
  range.mismatch.pheno.late=range.mismatch.pheno.late,
  ms.pred=length(range.mismatch.yr),
  ms.pred.pheno=length(range.mismatch.pheno.med))

#Parameters to track
my.parameters.joint<- c(##### AGE MODEL #####
                        #'age.at.B.est','mismatch.all', #removed because using too much memory
                        'int_age','beta_prim9','sd.prim','sd.year.prim',
                        ##### SURV MODEL #####
                        "phi.A", 'sigma.phiA', 'eps.phiA.t',
                        "phi.J", 'sigma.phiJ', 'eps.phiJ.t',
                        "int.r.ad", 'sigma.r.ad', 'eps.r.t.ad',
                        "int.r.juv", 'sigma.r.juv', 'eps.r.t.juv',
                        'int.p', 'sigma.p', 'eps.p.t',
                        "beta.p.L","beta.p.f","beta.ms","beta.peakN",
                        'p.low.m','p.low.f',
                        'p.high.m','p.high.f', 
                        'paH','paE',
                        'lpjH', 'lpjE', 'pjH', 'pjE', 'pjL',
                        'Sj.t',#'Sj.t.new',#'numbers.t.old','numbers.t.new',
                        'Sa', 'Sa.mean.p1', 'Sa.mean.p2', 'Sa.mean.p3', 'Sj.mean.p1', 'Sj.mean.p2', 'Sj.mean.p3', 
                        'mean_ms_all','sd_mean_ms_all',
                        'beta.hunt.ad','beta.hunt.juv',
                        #'Sj', #Don't forget to remove this if using individual-level covariate, this will take too much memory
                        #'mean.yr.ms','sd.yr.ms',
                        'r.ad', 'r.juv','S.pred.ms','S.pred.ms.bin','S.pred.ms.early','S.pred.ms.med','S.pred.ms.late')

my.parameters.joint.thin2<- c('logProb_age.at.B','logProb_y')

#Initial values
the.inits.joint <- function(){list(
  #### AGE MODEL ##### #
  int_age = 33,
  beta_prim9 = 3.15,
  sd.prim = 1.5,
  sd.year.prim = 1.5,
  #age.at.B.all = init.age.at.B.all,
  age.at.B.est = rep(30, times=length(b.age.unk.vec)),
  #gi = rep(c(25,(length(index.ms.values)+1)), times=c(length(juv),length(ad))),
  #age.at.B.est.bis = rep(30, times=length(ID.ms.unk)),
  
  ##### SURVIVAL MODEL ##### #
  phi.A= runif(1, -3, 3),
  phi.J= runif(1, -3, 3),
  beta.p.L = runif(1, 0, 1),
  beta.p.f = runif(1,-2, 2),
  beta.ms = runif(1,-2,2),
  beta.peakN = runif(1,-2,2),
  beta.hunt.ad = runif(my.constants.joint$n.phunt,-2,2),
  beta.hunt.juv = runif(my.constants.joint$n.phunt,-2,2),
  int.p = runif(1,-5,5),
  int.r.ad = runif(1, -5, 5),
  int.r.juv = runif(1, -5, 5),
  sigma.phiA = runif(1,0,3),
  sigma.phiJ = runif(1,0,3),
  sigma.p = runif(1,0,3),
  sigma.r.ad = runif(1,0,3),
  sigma.r.juv = runif(1,0,3),
  paH = runif(2,0,1),
  paE = runif(2,0,1),
  lpjH = runif(2,-5,5),
  lpjE = runif(2,-5,5))}


#Model
joint_mismatch_surv_analysis_optz<-nimbleCode({
  
  #################### AGE MODEL ######################
  
  #prior on beta 9th primary to determine age
  beta_prim9~dunif(0,50)
  int_age~dunif(0,50)
  sd.prim~dunif(0,5)
  sd.year.prim~dunif(0,5)
  
  #Random year intercepts for growth 'curve'
  for(t in 2:n.yrs){
    eps.ageB.t[t] ~ dnorm(mean=0, sd=sd.year.prim)
  }
  
  #Because no webtags put in 1990, impossible to estimate random intercept of that year for effect of 9prim on age, put 0 for now, not a big difference anyway
  eps.ageB.t[1]<-0
  
  for(i in 1:n.wbtg){
    
    # determine relationship between age and 9prim lgth (from known age juveniles)
    mu[i] <- int_age + beta_prim9 * prim9[b.age.kwn.vec[i]] + eps.ageB.t[yr.kwn[i]]
    
    age.at.B[i] ~ dnorm(mu[i], sd=sd.prim) #Likelihood
        
    ms0[b.age.kwn.vec[i]]<- mismatch.known[i] - min.ms + 1# calculate mismatch between known age and peak nitrogen in year of marking
    
  }
  
  #Predict age of juveniles of unknown age from modeled relationshio between age and 9th primary feather
  for(i in 1:n.not.wbtg){
    
    age.at.B.est[i] ~ dnorm(mu.pred[i], sd=sd.prim) # draw age of bird from normal distribution around linear predictor mean 
    
    mu.pred[i]<- int_age + beta_prim9 * prim9[b.age.unk.vec[i]] + eps.ageB.t[yr.unk[i]]# Linear predictor for age of bird i based on its prim9 lgth
    
    #age.at.B.all[b.age.unk.vec[i]] <- age.at.B.est[i] # Put age of bird in vector of ages
    
    ms0[b.age.unk.vec[i]] <- round((julian_date_B_unk[i]-age.at.B.est[i]) - peakN[yr.unk[i]]) - min.ms + 1# calculate mismatch between predicted age and peak nitrogen in year of marking, and put at right place in vector of mismatch values
    
  }
  
  for(i in 1:n.juv){
    ms.class[i]<-ms.class.vec[ms0[i]]  # Classify mismatch in reduced pre-determined bins
  }
  
  #Keep a sample of estimated ages to see how well they are estimated but not the whole dataset because it takes too much memory
  # for(i in 1:n.age.save){
  #   
  #   ms0.keep[i]<-ms0[age.keep[i]]
  #   
  # }
  
  ############################################################################################# #
  ########################################## SURVIVAL MODEL  #################################### 
  ############################################################################################# #
  
  # ------------------------------------------------
  # Main parameters:
  # Sj: juvenile survival probability
  # Sa: adult survival probability
  
  # r.ad: recovery probability for adults
  # r.juv: recovery probability for juveniles
  # p.high.m: recapture probability for males in highly-observable group
  # p.high.f: recapture probability for females in highly-observable group
  # p.low.m: recapture probability for males in low observable group
  # p.low.f: recapture probability for females in low observable group

  # paH: Probability for adults to be in high encounter group (for initial states matrix)
  # paE: Probability for adults to emigrate permanently from study area (for transition matrix)
  # pjH: Probability for juveniles to be part of high encounter group conditional on 1st year survival (for transition matrix)
  # pjE: Probability for juveniles to emigrate permanently from study area conditional on 1st year survival (for transition matrix)
  
  # ------------------------------------------------
  # States (S):
  # 1 alive young
  # 2 alive adult and highly observable
  # 3 alive adult and weakly observable
  # 4 alive adult and permanently emigrated
  # 5 Newly dead and recovered
  # 6 Recently dead, but not recovered, or dead (absorbing)
  
  # Observations (O):
  # 1 seen alive
  # 2 recovered dead (coded in transition matrix following Kéry & Schaub 2010)
  # 3 neither seen nor recovered
  # ------------------------------------------------
  
  
  ########################################################## #
  ###################  CONSTRAINTS  ########################
  ########################################################## #
  
  #####Survival adults#### #
  
  #First period: no beta.hunt
  for (t in 1:end.occasion.p1){
    #Constraints for survival: estimate yearly survival rate
    logit(Sa[t]) <- phi.A + eps.phiA.t[t] #Adults
  }
  
  #2nd (and/or 3rd) periods: go from last occasion of first period +1 until end of dataset
  for (t in (end.occasion.p1+1):(n.occasions-1)){#On part de la première occasion après le début de la chasse et on va à la fin des données
    #Constraints for survival: estimate yearly survival rate
    logit(Sa[t]) <- phi.A + eps.phiA.t[t] + beta.hunt.ad[hunt[t]] #Adults
  }
  
  #####Survival juveniles#### #
  
  
  #Index juvenile survival by mismatch class rather than by individual (results in a lot less transition matrix slices, and saves a lot of memory and model graph building time - days)
  for(c in 1:n.ms.class){
    
    for (t in 1:end.occasion.p1){# période sans chasse
      logit(Sj[c,t])<- phi.J + beta.ms * mismatch.scaled[c] + beta.peakN * scaled.peakN[t] + eps.phiJ.t[t]
    }#for (t in 1:end.occasion.p1)
    
    for (t in (end.occasion.p1+1):(n.occasions-1)){#période avec chasse
      logit(Sj[c,t])<- phi.J + beta.ms * mismatch.scaled[c] + beta.hunt.juv[hunt[t]] + beta.peakN * scaled.peakN[t] + eps.phiJ.t[t]
    }#for (t in (end.occasion.p1+1):(n.occasions-1))
  }#for(c in 1:n.ms.class){
  
  #### EVENTS #### #
  for (t in 2:n.occasions){
    
    #Capture probability (here variable through time (random effect))
    logit(p.high.m[t]) <- int.p + eps.p.t[t]              # P.cap for males
    logit(p.high.f[t]) <- int.p + eps.p.t[t] + beta.p.f   # P.cap for females
    
    #New parametrization to constrain p.low<p.high
    p.low.m[t] <- p.high.m[t] * beta.p.L 
    p.low.f[t] <- p.high.f[t] * beta.p.L
    
    #Recovery probability (here constrained equal through time)
    logit(r.ad[t])  <- int.r.ad  + eps.r.t.ad[t]  
    logit(r.juv[t])  <- int.r.juv  + eps.r.t.juv[t]  
  }#for (t in 2:n.occasions)
  
  ########################################################## #
  #####################  PRIORS  ###########################
  ########################################################## #
  
  
  ## heterogeneity parameters ##
  for(s in 1:n.sex){
    paH[s] ~ dunif(0,1) #Probability of being in high observability group for adults
    paE[s] ~ dunif(0,1) #Probability of emigrating for adults
    
    #Juveniles 
    #Loop twice for males and females
    #Males
    lpjH[s] ~ dnorm(0, sd = 1.6) #High-p
    lpjE[s] ~ dnorm(0, sd = 1.6) #Emigration
    
    # constrain the transitions such that the sum of probabilities for surviving juveniles going to adult high, low and emigrated is = 1
    pjH[s] <- exp(lpjH[s]) / (1 + exp(lpjH[s]) + exp(lpjE[s])) #High-p
    pjE[s] <- exp(lpjE[s]) / (1 + exp(lpjH[s]) + exp(lpjE[s])) #Emigration
    
    # last transition probability
    pjL[s] <- 1 - pjH[s] - pjE[s]                          #Low-p, constrain so that all 3 add up to 1
  }#Fin for s in 1:sex
  
  #Betas
  beta.p.L ~ dunif(0,1) #Prior for mean effect of heterogeneity on p, constrained to be function of p.high, so that p.low<p.high (i.e., no jumping between both parameters during sampling)
  beta.p.f ~ dnorm(0,sd=1.6) #Prior for effect of female on recap, should be 0 if heterogeneity is well set up
  beta.peakN ~ dnorm (0,sd=1.6) #Prior for effect of annual peak Nitrogen date on juvenile survival
  
  #mismatch effect
  beta.ms ~ dnorm(0, sd=1.6)
  
  #hunt effect
  for(ph in 1:n.phunt){
    beta.hunt.ad[ph] ~ dnorm(0, sd=1.6)
    beta.hunt.juv[ph] ~ dnorm(0, sd=1.6)
  }#for(ph in 1:n.phunt)
  
  #Random effects
  phi.A ~ dnorm(0, sd=1.6) # Intercept survival adults
  phi.J ~ dnorm(0, sd=1.6) # Intercept survival juveniles
  
  sigma.phiA ~ dunif(0,5)  # Variance survival adults
  sigma.phiJ ~ dunif(0,5)  # Variance survival juveniles
  
  int.p ~ dnorm (0, sd=1.6) # Intercept recap
  sigma.p ~ dunif(0,5)     # Variance recap
  
  int.r.ad ~ dnorm(0, sd=1.6)  # Intercept recoveries adults
  sigma.r.ad ~ dunif(0,5)     # Variance recoveries adults
  
  int.r.juv ~ dnorm(0, sd=1.6)  # Intercept recoveries juveniles
  sigma.r.juv ~ dunif(0,5)     # Variance recoveries juveniles
  
  # Random time intercepts
  for(t in 2:n.occasions){
    eps.phiA.t[t-1] ~ dnorm (0, sd=sigma.phiA)    # Prior survival adult random intercept
    
    eps.phiJ.t[t-1] ~ dnorm (0, sd=sigma.phiJ)    # Prior survival juvenile random intercept
    
    eps.p.t[t] ~ dnorm (0, sd=sigma.p)         # Prior recapture random intercept
    
    eps.r.t.ad[t] ~ dnorm (0, sd=sigma.r.ad)         # Prior recovery adults random intercept
    
    eps.r.t.juv[t] ~ dnorm (0, sd=sigma.r.juv)         # Prior recovery juveniles random intercept
    
  }#for(t in 2:n.occasions)
  
  
  ########################################################## #
  #####################  MATRICES  #########################
  ########################################################## #
  
  ###################  Initial states  ######################
  for(s in 1:n.sex){
  delta[1,1,s] <- 1            #Probability of being young upon first capture = 1 for i in juv
  delta[1,2,s] <- 0            #Probability of being adult high-p upon first capture = 0
  delta[1,3,s] <- 0            #Probability of being adult low-p upon first capture = 0
  delta[1,4,s] <- 0            #Probability of being adult emigrated upon first capture = 0
  delta[1,5,s] <- 0
  delta[1,6,s] <- 0
  
  
  
  
  #Adults
  delta[2,1,s] <- 0           #Probability of being young upon first capture = 0 for i in ad
  delta[2,2,s] <- paH[s]        #Probability of being adult high-p upon first capture
  delta[2,3,s] <- 1-paH[s]      #Probability of being adult low-p upon first capture
  delta[2,4,s] <- 0           #Probability of being adult emigrated upon first capture = 0
  delta[2,5,s] <- 0
  delta[2,6,s] <- 0
  }
  
  
  
  
  
  ###################  TRANSITIONS  #########################
  
  for(c in 1:n.ms.class){ #Indexing by classes of mismatch rather than by individual to save memory and computing time
    for(s in 1:n.sex){
      for(t in 1:n.occasions-1){
        #Alive young
        gamma[1,1,c,t,s] <- 0                           #Never stays alive young
        gamma[1,2,c,t,s] <- Sj[c,t] * pjH[s]               #Can survive and become high-p adult
        gamma[1,3,c,t,s] <- Sj[c,t] * pjL[s]               #Can survive and become low-p adult
        gamma[1,4,c,t,s] <- Sj[c,t] * pjE[s]               #Can survive and emigrate
        gamma[1,5,c,t,s] <- (1-Sj[c,t]) *r.juv[t+1]         #Can die and be recovered
        gamma[1,6,c,t,s] <- (1-Sj[c,t]) *(1-r.juv[t+1])     #Can die and not be recovered
        
        #Alive adults (high-p)
        gamma[2,1,c,t,s] <- 0                     #Never becomes alive young
        gamma[2,2,c,t,s] <- Sa[t] * (1-paE[s])    #Survival (stays high-p)
        gamma[2,3,c,t,s] <- 0                     #Survival (no tr. to low-p)
        gamma[2,4,c,t,s] <- Sa[t] * paE[s]        #Survival (emigrates)
        gamma[2,5,c,t,s] <- (1-Sa[t]) *r.ad[t+1]     #Can die and be recovered
        gamma[2,6,c,t,s] <- (1-Sa[t]) *(1-r.ad[t+1]) #Can die and not be recovered
        
        #Alive adults (low-p)
        gamma[3,1,c,t,s] <- 0                     #Never becomes alive young
        gamma[3,2,c,t,s] <- 0                     #Survival (no tr. to high-p)
        gamma[3,3,c,t,s] <- Sa[t] * (1-paE[s])    #Survival (stays low-p)
        gamma[3,4,c,t,s] <- Sa[t] * paE[s]        #Survival (emigrates)
        gamma[3,5,c,t,s] <- (1-Sa[t]) *r.ad[t+1]     #Can die and be recovered
        gamma[3,6,c,t,s] <- (1-Sa[t]) *(1-r.ad[t+1]) #Can die and not be recovered
        
        #Alive adults (low-p)
        gamma[4,1,c,t,s] <- 0                     #Never becomes alive young
        gamma[4,2,c,t,s] <- 0                     #Survival (no tr. to high-p)
        gamma[4,3,c,t,s] <- 0                     #Survival (no tr. to low-p)
        gamma[4,4,c,t,s] <- Sa[t]                 #Survival (stays emigrated)
        gamma[4,5,c,t,s] <- (1-Sa[t]) *r.ad[t+1]     #Can die and be recovered
        gamma[4,6,c,t,s] <- (1-Sa[t]) *(1-r.ad[t+1]) #Can die and not be recovered
        
        #newly dead (all ages)
        gamma[5,1,c,t,s] <- 0
        gamma[5,2,c,t,s] <- 0
        gamma[5,3,c,t,s] <- 0
        gamma[5,4,c,t,s] <- 0
        gamma[5,5,c,t,s] <- 0
        gamma[5,6,c,t,s] <- 1
        
        #"Old" dead (all ages) (or newly dead and not recovered)
        gamma[6,1,c,t,s] <- 0
        gamma[6,2,c,t,s] <- 0
        gamma[6,3,c,t,s] <- 0
        gamma[6,4,c,t,s] <- 0
        gamma[6,5,c,t,s] <- 0
        gamma[6,6,c,t,s] <- 1
        
      }# for t in n.occasions-1
    }#for s in 1:n.sex
  }# for c in ms.classes
  
  
  
  ##################  SURVIVAL  #######################
  
  #for adults, only need to index gamma by time and not individual since they share the same matrix slice for a given time step (no individual covariates for AD, saves on model graph building time)
  #gi.ad takes the value of number of juveniles +1,
  for(s in 1:n.sex){
    for (t in 1:(n.occasions-1)){
      
      #Alive young
      gamma[1,1,gi.ad,t,s] <- 0                      #Never stays alive young
      gamma[1,2,gi.ad,t,s] <- 0                      #No juvenile survival for adult birds  
      gamma[1,3,gi.ad,t,s] <- 0                      #No juvenile survival for adult birds  
      gamma[1,4,gi.ad,t,s] <- 0                      #No juvenile survival for adult birds  
      gamma[1,5,gi.ad,t,s] <- 0                      #No juvenile survival for adult birds  
      gamma[1,6,gi.ad,t,s] <- 0                      #No juvenile survival for adult birds  
      
      #Alive adults (high-p)
      gamma[2,1,gi.ad,t,s] <- 0                     #Never becomes alive young
      gamma[2,2,gi.ad,t,s] <- Sa[t] * (1-paE[s])    #Survival (stays high-p) 
      gamma[2,3,gi.ad,t,s] <- 0                     #Survival (no tr. to low-p)
      gamma[2,4,gi.ad,t,s] <- Sa[t] * paE[s]        #Survival (emigrates)
      gamma[2,5,gi.ad,t,s] <- (1-Sa[t]) *r.ad[t+1]     #Can die and be recovered 
      gamma[2,6,gi.ad,t,s] <- (1-Sa[t]) *(1-r.ad[t+1]) #Can die and not be recovered 
      
      #Alive adults (low-p)
      gamma[3,1,gi.ad,t,s] <- 0                     #Never becomes alive young
      gamma[3,2,gi.ad,t,s] <- 0                     #Survival (no tr. to high-p) 
      gamma[3,3,gi.ad,t,s] <- Sa[t] * (1-paE[s])    #Survival (stays low-p) 
      gamma[3,4,gi.ad,t,s] <- Sa[t] * paE[s]        #Survival (emigrates)
      gamma[3,5,gi.ad,t,s] <- (1-Sa[t]) *r.ad[t+1]     #Can die and be recovered 
      gamma[3,6,gi.ad,t,s] <- (1-Sa[t]) *(1-r.ad[t+1]) #Can die and not be recovered 
      
      #Alive adults (low-p)
      gamma[4,1,gi.ad,t,s] <- 0                     #Never becomes alive young
      gamma[4,2,gi.ad,t,s] <- 0                     #Survival (no tr. to high-p) 
      gamma[4,3,gi.ad,t,s] <- 0                     #Survival (no tr. to low-p) 
      gamma[4,4,gi.ad,t,s] <- Sa[t]                 #Survival (stays emigrated) 
      gamma[4,5,gi.ad,t,s] <- (1-Sa[t]) *r.ad[t+1]     #Can die and be recovered 
      gamma[4,6,gi.ad,t,s] <- (1-Sa[t]) *(1-r.ad[t+1]) #Can die and not be recovered 
      
      #newly dead (all ages)
      gamma[5,1,gi.ad,t,s] <- 0
      gamma[5,2,gi.ad,t,s] <- 0
      gamma[5,3,gi.ad,t,s] <- 0
      gamma[5,4,gi.ad,t,s] <- 0
      gamma[5,5,gi.ad,t,s] <- 0
      gamma[5,6,gi.ad,t,s] <- 1
      
      #"Old" dead (all ages) (or newly dead and not recovered)
      gamma[6,1,gi.ad,t,s] <- 0
      gamma[6,2,gi.ad,t,s] <- 0
      gamma[6,3,gi.ad,t,s] <- 0
      gamma[6,4,gi.ad,t,s] <- 0
      gamma[6,5,gi.ad,t,s] <- 0
      gamma[6,6,gi.ad,t,s] <- 1
      
    }# for t in f[i]:(n.occasions-1)
  }#for(s in 1:n.sex)
  ##################  EVENTS  ####################
  #Events with capture probability for males first and females 2nd
  for(m in 1:n.occasions){ # m takes the value of the marking occasion, because all individuals marked in the same occasion share the same capture and recovery probabilities
    
    for (t in (m+1):n.occasions){
      # Define probabilities of O(t) given S(t)
      
      #Alive young never observed as young beyond 1st capture
      omega[1,1,m,t,1] <- 0  
      omega[1,2,m,t,1] <- 0
      omega[1,3,m,t,1] <- 1  
      
      #Alive adult low.p can be observed 
      omega[2,1,m,t,1] <- p.high.m[t]
      omega[2,2,m,t,1] <- 0
      omega[2,3,m,t,1] <- 1 - p.high.m[t]
      
      #Alive and adult can be observed (to be split in two high and low observability classes)
      omega[3,1,m,t,1] <- p.low.m[t]
      omega[3,2,m,t,1] <- 0
      omega[3,3,m,t,1] <- 1 - p.low.m[t]
      
      #Alive and adult emigrated p = 0
      omega[4,1,m,t,1] <- 0
      omega[4,2,m,t,1] <- 0
      omega[4,3,m,t,1] <- 1
      
      omega[5,1,m,t,1] <- 0
      omega[5,2,m,t,1] <- 1 #recoveries coded in transitions so here this = 1
      omega[5,3,m,t,1] <- 0
      
      omega[6,1,m,t,1] <- 0
      omega[6,2,m,t,1] <- 0
      omega[6,3,m,t,1] <- 1
      
      
      
    } # for t in (f[i]+1):(n.occasions)
    
    #Here, omega is indexed by marking occasion (3rd dim) and we need to give values for initial capture (at time t = m)
    
    #Alive young observed as young first capture
    omega[1,1,m,m,1] <- 1  
    omega[1,2,m,m,1] <- 0
    omega[1,3,m,m,1] <- 0  
    
    #Alive and adult observed on first capture
    omega[2,1,m,m,1] <- 1
    omega[2,2,m,m,1] <- 0
    omega[2,3,m,m,1] <- 0
    
    #Alive and adult observed on first capture
    omega[3,1,m,m,1] <- 1
    omega[3,2,m,m,1] <- 0
    omega[3,3,m,m,1] <- 0
    
    #Alive and emigrated impossible on first capture
    omega[4,1,m,m,1] <- 1
    omega[4,2,m,m,1] <- 0
    omega[4,3,m,m,1] <- 0
    
    #Newly dead impossible on first capture
    omega[5,1,m,m,1] <- 0
    omega[5,2,m,m,1] <- 0 
    omega[5,3,m,m,1] <- 0 
    
    #Dead impossible on first capture
    omega[6,1,m,m,1] <- 0
    omega[6,2,m,m,1] <- 0 
    omega[6,3,m,m,1] <- 0 
    
    
  }#for m in 1:n.occasions
  
  #Events with capture probability for females 2nd (last array dimension is sex)
  for(m in 1:n.occasions){ # m takes the value of the marking occasion, because all individuals marked in the same occasion share the same capture an recovery probabilities
    
    for (t in (m+1):n.occasions){
      # Define probabilities of O(t) given S(t)
      
      #Alive young never observed as young beyond 1st capture
      omega[1,1,m,t,2] <- 0  
      omega[1,2,m,t,2] <- 0
      omega[1,3,m,t,2] <- 1  
      
      #Alive adult low.p can be observed 
      omega[2,1,m,t,2] <- p.high.f[t]
      omega[2,2,m,t,2] <- 0
      omega[2,3,m,t,2] <- 1 - p.high.f[t]
      
      #Alive and adult can be observed (to be split in two high and low observability classes)
      omega[3,1,m,t,2] <- p.low.f[t]
      omega[3,2,m,t,2] <- 0
      omega[3,3,m,t,2] <- 1 - p.low.f[t]
      
      #Alive and adult emigrated p = 0
      omega[4,1,m,t,2] <- 0
      omega[4,2,m,t,2] <- 0
      omega[4,3,m,t,2] <- 1
      
      omega[5,1,m,t,2] <- 0
      omega[5,2,m,t,2] <- 1 #recoveries coded in transitions so here this = 1
      omega[5,3,m,t,2] <- 0
      
      omega[6,1,m,t,2] <- 0
      omega[6,2,m,t,2] <- 0
      omega[6,3,m,t,2] <- 1
      
      
      
    } # for t in (f[i]+1):(n.occasions)
    
    #Here, omega is indexed by marking occasion (3rd dim) and we need to give values for initial capture (at time = m)
    
    #Alive young observed as young first capture
    omega[1,1,m,m,2] <- 1  
    omega[1,2,m,m,2] <- 0
    omega[1,3,m,m,2] <- 0  
    
    #Alive and adult observed on first capture
    omega[2,1,m,m,2] <- 1
    omega[2,2,m,m,2] <- 0
    omega[2,3,m,m,2] <- 0
    
    #Alive and adult observed on first capture
    omega[3,1,m,m,2] <- 1
    omega[3,2,m,m,2] <- 0
    omega[3,3,m,m,2] <- 0
    
    #Alive and emigrated impossible on first capture
    omega[4,1,m,m,2] <- 0
    omega[4,2,m,m,2] <- 0
    omega[4,3,m,m,2] <- 0
    
    #Newly dead impossible on first capture
    omega[5,1,m,m,2] <- 0
    omega[5,2,m,m,2] <- 0 
    omega[5,3,m,m,2] <- 0 
    
    #Dead impossible on first capture
    omega[6,1,m,m,2] <- 0
    omega[6,2,m,m,2] <- 0 
    omega[6,3,m,m,2] <- 0 
    
    
  }#for m in 1:n.occasions

  

  ########################################################## #
  ###################### LIKELIHOOD ########################
  ########################################################## #

  # To associate correct transition matrix slice to adults in likelihood statement below
  ms.class[(n.juv+1):nind]<-n.ms.class+1
  
  # Likelihood
  for (i in 1:nind) {
    # note: using (f[i] + 1):n.occasions below generates errors
    y[i, f[i]:n.occasions] ~ dDHMMo(init = delta[age[i],1:6,sex[i]], 
                                    probTrans = gamma[1:6,1:6,ms.class[i],f[i]:(n.occasions-1),sex[i]],
                                    probObs = omega[1:6,1:3,f[i],f[i]:n.occasions,sex[i]],
                                    len = length(f[i]:n.occasions),
                                    checkRowSums = 0) 
    
  }
  
  #Derived quantity: average survival
  #Derived quantity: average survival
  Sa.mean.p1 <- mean(Sa[1:end.occasion.p1])
  Sj.mean.p1 <- mean(Sj.t[2:end.occasion.p1])
  
  Sa.mean.p2 <- mean(Sa[(end.occasion.p1+1):end.occasion.p2])
  Sj.mean.p2 <- mean(Sj.t[(end.occasion.p1+1):end.occasion.p2])
  
  Sa.mean.p3 <- mean(Sa[(end.occasion.p2+1):(n.occasions-1)])
  Sj.mean.p3 <- mean(Sj.t[(end.occasion.p2+1):(n.occasions-1)])
  
  
  #Calculate survival for each year based on average mismatch, with only known mismatch values
  # for (t in 1:end.occasion.p1){# période sans chasse
  #   
  #   mean_ms[t]<-(average_yearly_mismatch[t]-scale.ms.mean)/scale.ms.sd
  #   
  #   logit(Sj.t[t])<-  phi.J + beta.ms * mean_ms[t] + eps.phiJ.t[t]
  #   
  # }
  # 
  # for (t in (end.occasion.p1+1):(n.occasions-1)){#périodes avec chasse
  #   
  #   mean_ms[t]<-(average_yearly_mismatch[t]-scale.ms.mean)/scale.ms.sd
  #   
  #   logit(Sj.t[t])<-  phi.J + beta.ms * mean_ms[t] + eps.phiJ.t[t] + beta.hunt.juv[hunt[t]]
  # }
  
  #Calculate survival for each year based on average mismatch, including estimated mismatch values
  for(t in 2:(n.occasions-1)){
  
  mean_ms_all[t]<- calc_avg_ms0_t(ms0_vec=ms0[1:n.juv], ms0_avg_index_vec=juv.f.indices[1:n.juv.marked.yr[t],t])
  
  #Derived parameter: sd of estimated and real mismatches to have an idea of the variability within years
  sd_mean_ms_all[t]<-calc_sd_ms0_t(ms0_vec=ms0[1:n.juv], ms0_avg_index_vec=juv.f.indices[1:n.juv.marked.yr[t],t])
  
  }
  
  for (t in 2:end.occasion.p1){# période sans chasse
    
    mean_ms_all_scaled[t]<-(mean_ms_all[t]+min.ms-1-scale.ms.mean)/scale.ms.sd
    
    logit(Sj.t[t])<-  phi.J + beta.ms * mean_ms_all_scaled[t] + beta.peakN * scaled.peakN[t] + eps.phiJ.t[t]
    
  }
  
  for (t in (end.occasion.p1+1):(n.occasions-1)){#périodes avec chasse
    
    mean_ms_all_scaled[t]<-(mean_ms_all[t]+min.ms-1-scale.ms.mean)/scale.ms.sd
    
    logit(Sj.t[t])<-  phi.J + beta.ms * mean_ms_all_scaled[t] +  beta.hunt.juv[hunt[t]] + beta.peakN * scaled.peakN[t] + eps.phiJ.t[t]
  }
  
  
  
  
  #mismatch effect on S - continuous predictions
  for(i in 1:ms.pred){
    
    logit(S.pred.ms[i])<- phi.J + beta.ms * range.mismatch[i]
    
  }
  
  #mismatch effect on S - predictions including peakN effect, with range of mismatch observed for each pheno bin, not really interesting in this model because no interaction.
  for(i in 1:ms.pred.pheno){
    
    logit(S.pred.ms.early[i])<- phi.J + beta.ms * range.mismatch.pheno.early[i] + beta.peakN * scaled.pN.early
    logit(S.pred.ms.med[i])<- phi.J + beta.ms * range.mismatch.pheno.med[i] + beta.peakN * scaled.pN.med
    logit(S.pred.ms.late[i])<- phi.J + beta.ms * range.mismatch.pheno.late[i] + beta.peakN * scaled.pN.late
    
  }
  
  
  #mismatch effect on S - prediction for each bin
  for(c in 1:n.ms.class){
    
    S.pred.ms.bin[c]<-phi.J + beta.ms * mismatch.scaled[c]
    
  }
})



# nimbleFunction to calculate annual average juvenile survival in model:
calc_avg_ms0_t <- nimbleFunction(
  run = function(ms0_vec = double(1), ms0_avg_index_vec = double(1)) { # all model variables should be doubles, even if they hold integer values
    ans <- 0
    n <- length(ms0_avg_index_vec)
    if(n==0) return(ans)
    for(i in 1:n) {
      ans <- ans + (ms0_vec[ms0_avg_index_vec[i]])
    }
    ans <- ans/n
    return(ans)
    returnType(double())
  })


calc_sd_ms0_t <- nimbleFunction(
  run = function(ms0_vec = double(1), ms0_avg_index_vec = double(1)) { # all model variables should be doubles, even if they hold integer values
    mn <- 0
    ans.sd <- 0
    n <- length(ms0_avg_index_vec)
    if(n==0) return(ans.sd)
    for(i in 1:n) {
      mn <- mn + (ms0_vec[ms0_avg_index_vec[i]])
    }
    
    mn<-mn/n
    
    SSE <- 0
    for(i in 1:n) {
      SSE <- SSE + ((ms0_vec[ms0_avg_index_vec[i]]) - mn)^2 
    }
    ans.sd <- (SSE/(n-1))^(1/2)
    return(ans.sd)
    returnType(double())
  })




myConstants = my.constants.joint
myData = my.data.joint
myParameters = my.parameters.joint
myParameters2 = my.parameters.joint.thin2
myCode = joint_mismatch_surv_analysis_optz
b.age.unk.vec=b.age.unk.vec




Sys.time()
#### CALCUL CANADA MULTIPLE SAVE CODE

debut=Sys.time()
# Build model
gsg.mod <- nimbleModel(code = myCode,
                       data = myData,
                       constants = myConstants,
                       inits = the.inits.joint(),
                       calculate=F)

paste('built model',Sys.time())
(time_built=Sys.time()-debut)

# build MCMC
debut=Sys.time()
conf <- configureMCMC(gsg.mod,monitors=myParameters, monitors2=myParameters2, 
                      useConjugacy = F, enableWAIC = T)
Sys.time()
gsg.MCMC <- buildMCMC(conf)
paste('built MCMC',Sys.time())

# Compile model and MCMC 
C.gsg.mod<-compileNimble(gsg.mod)
Cgsg.MCMC <- compileNimble(gsg.MCMC, project=gsg.mod)
paste('compiled MCMC and model',Sys.time())
fin=Sys.time()
(time_mcmc=fin-debut)


# Run model in many smaller runs to save progress as we go and avoid losing information in case of cluster shutdown (yes, this happened)

start.chain=Sys.time()
#First run, keep all samples for the first tests, can always discard them later
Cgsg.MCMC$run(niter=30, nburnin=n.burnin, thin=n.thin, thin2=n.thin2,  progressBar=T)
end.chain=Sys.time()
(time.chain<-end.chain-start.chain)

#Extract and save samples
chain1<-as.matrix(Cgsg.MCMC$mvSamples)
saveRDS(chain1,file = './result_chain1/chain1_run0.rds')

cat(paste0('saved_samples init, run ',i,' at ',Sys.time(),'\n', sep=''))

#Save log probabilities to compute WAIC afterwards, and be able to remove parameters from age model (some pwaic>0.4 in age model, which make the WAIC unstable and unreliable)
#Should be OK since only interested in comparing different survival models and age model never changes, so no need to include its parameters for model selection.
chain1_logProbs<-as.matrix(Cgsg.MCMC$mvSamples2)
cat(paste0('matrix logprobs init, run ',i,' at ',Sys.time(),'\n', sep=''))

saveRDS(chain1_logProbs,file = './result_chain1/chain1_run0_logProbs.rds')
cat(paste0('saved logprobs init, run ',i,' at ',Sys.time(),'\n', sep=''))


#Extract model and MCMC state and save them
stateList <- list(modelState = getModelState(C.gsg.mod),
                  mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                  rs = .Random.seed)

saveRDS(stateList, file = './result_chain1/model_state.rds')

#Extract WAIC up to here
WAIC<-Cgsg.MCMC$getWAIC()

rm(stateList,WAIC)

n.runs<-9

for(i in 1:n.runs){
start.chain=Sys.time()
#################################################################################################################  
#Run 1.1
Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
end.chain=Sys.time()
(time.chain<-end.chain-start.chain)

cat(paste0('ran chain 1/5, run ',i,' at ',Sys.time(),'\n', sep=''))

#Extract and save samples
chain1<-as.matrix(Cgsg.MCMC$mvSamples)
saveRDS(chain1,file = paste0('./result_chain1/chain1_run',i,'.rds', sep=''))
cat(paste0('saved_samples 1/5, run ',i,' at ',Sys.time(),'\n', sep=''))

chain1_logProbs<-as.matrix(Cgsg.MCMC$mvSamples2)
cat(paste0('matrix logProbs 1/5, run ',i,' at ',Sys.time(),'\n', sep=''))

saveRDS(chain1_logProbs,file = paste0('./result_chain1/chain1_run',i,'_logProbs.rds', sep=''))
cat(paste0('saved logProbs 1/5, run ',i,' at ',Sys.time(),'\n', sep=''))


#Extract model and MCMC state and save them
stateList <- list(modelState = getModelState(C.gsg.mod),
                  mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                  rs = .Random.seed)

saveRDS(stateList, file = './result_chain1/model_state.rds')
cat(paste0('saved statelist 1/5, run ',i,' at ',Sys.time(),'\n', sep=''))

#Extract WAIC up to here
WAIC<-Cgsg.MCMC$getWAIC()

saveRDS(WAIC, file=paste0('./result_chain1/waic_run',i,'.rds', sep=''))
cat(paste0('saved WAIC 1/5, run ',i,' at ',Sys.time(),'\n', sep=''))

rm(stateList,WAIC)

#################################################################################################################  
#Run 1.2
Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)

cat(paste0('ran chain 2/5, run ',i,' at ',Sys.time(),'\n', sep=''))


#Extract and save samples
chain1<-rbind(chain1,as.matrix(Cgsg.MCMC$mvSamples))
saveRDS(chain1,file = paste0('./result_chain1/chain1_run',i,'.rds', sep=''))
cat(paste0('saved_samples 2/5, run ',i,' at ',Sys.time(),'\n', sep=''))

chain1_logProbs<-rbind(chain1_logProbs,as.matrix(Cgsg.MCMC$mvSamples2))
cat(paste0('matrix logProbs 2/5, run ',i,' at ',Sys.time(),'\n', sep=''))

saveRDS(chain1_logProbs,file = paste0('./result_chain1/chain1_run',i,'_logProbs.rds', sep=''))
cat(paste0('saved logProbs 2/5, run ',i,' at ',Sys.time(),'\n', sep=''))


#Extract model and MCMC state and save them
stateList <- list(modelState = getModelState(C.gsg.mod),
                  mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                  rs = .Random.seed)

saveRDS(stateList, file = './result_chain1/model_state.rds')
cat(paste0('saved statelist 2/5, run ',i,' at ',Sys.time(),'\n', sep=''))

#Extract WAIC up to here
WAIC<-Cgsg.MCMC$getWAIC()

saveRDS(WAIC, file=paste0('./result_chain1/waic_run',i,'.rds', sep=''))
cat(paste0('saved WAIC 2/5, run ',i,' at ',Sys.time(),'\n', sep=''))

rm(stateList, WAIC)

#################################################################################################################  
#Run 1.3
Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
cat(paste0('ran chain 3/5, run ',i,' at ',Sys.time(),'\n', sep=''))

#Extract and save samples
chain1<-rbind(chain1,as.matrix(Cgsg.MCMC$mvSamples))
saveRDS(chain1,file = paste0('./result_chain1/chain1_run',i,'.rds', sep=''))
cat(paste0('saved_samples 3/5, run ',i,' at ',Sys.time(),'\n', sep=''))

chain1_logProbs<-rbind(chain1_logProbs,as.matrix(Cgsg.MCMC$mvSamples2))
cat(paste0('matrix logProbs 3/5, run ',i,' at ',Sys.time(),'\n', sep=''))

saveRDS(chain1_logProbs,file = paste0('./result_chain1/chain1_run',i,'_logProbs.rds', sep=''))
cat(paste0('saved logProbs 3/5, run ',i,' at ',Sys.time(),'\n', sep=''))


#Extract model and MCMC state and save them
stateList <- list(modelState = getModelState(C.gsg.mod),
                  mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                  rs = .Random.seed)

saveRDS(stateList, file = './result_chain1/model_state.rds')
cat(paste0('saved statelist 3/5, run ',i,' at ',Sys.time(),'\n', sep=''))

#Extract WAIC up to here
WAIC<-Cgsg.MCMC$getWAIC()

saveRDS(WAIC, file=paste0('./result_chain1/waic_run',i,'.rds', sep=''))
cat(paste0('saved WAIC 3/5, run ',i,' at ',Sys.time(),'\n', sep=''))

rm(stateList, WAIC)

#################################################################################################################  
#Run 1.4
Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
cat(paste0('ran chain 4/5, run ',i,' at ',Sys.time(),'\n', sep=''))

#Extract and save samples
chain1<-rbind(chain1,as.matrix(Cgsg.MCMC$mvSamples))
saveRDS(chain1,file = paste0('./result_chain1/chain1_run',i,'.rds', sep=''))
cat(paste0('saved_samples 4/5, run ',i,' at ',Sys.time(),'\n', sep=''))

chain1_logProbs<-rbind(chain1_logProbs,as.matrix(Cgsg.MCMC$mvSamples2))
cat(paste0('matrix logProbs 4/5, run ',i,' at ',Sys.time(),'\n', sep=''))

saveRDS(chain1_logProbs,file = paste0('./result_chain1/chain1_run',i,'_logProbs.rds', sep=''))
cat(paste0('saved logProbs 4/5, run ',i,' at ',Sys.time(),'\n', sep=''))


#Extract model and MCMC state and save them
stateList <- list(modelState = getModelState(C.gsg.mod),
                  mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                  rs = .Random.seed)

saveRDS(stateList, file = './result_chain1/model_state.rds')
cat(paste0('saved statelist 4/5, run ',i,' at ',Sys.time(),'\n', sep=''))

#Extract WAIC up to here
WAIC<-Cgsg.MCMC$getWAIC()

saveRDS(WAIC, file=paste0('./result_chain1/waic_run',i,'.rds', sep=''))
cat(paste0('saved WAIC 4/5, run ',i,' at ',Sys.time(),'\n', sep=''))

rm(stateList, WAIC)

#################################################################################################################  
#Run 1.5
Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
cat(paste0('ran chain 5/5, run ',i,' at ',Sys.time(),'\n', sep=''))

#Extract and save samples
chain1<-rbind(chain1,as.matrix(Cgsg.MCMC$mvSamples))
saveRDS(chain1,file = paste0('./result_chain1/chain1_run',i,'.rds', sep=''))
cat(paste0('saved_samples 5/5, run ',i,' at ',Sys.time(),'\n', sep=''))

chain1_logProbs<-rbind(chain1_logProbs,as.matrix(Cgsg.MCMC$mvSamples2))
cat(paste0('matrix logProbs 5/5, run ',i,' at ',Sys.time(),'\n', sep=''))

saveRDS(chain1_logProbs,file = paste0('./result_chain1/chain1_run',i,'_logProbs.rds', sep=''))
cat(paste0('saved logProbs 5/5, run ',i,' at ',Sys.time(),'\n', sep=''))


#Extract model and MCMC state and save them
stateList <- list(modelState = getModelState(C.gsg.mod),
                  mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                  rs = .Random.seed)

saveRDS(stateList, file = './result_chain1/model_state.rds')
cat(paste0('saved statelist 5/5, run ',i,' at ',Sys.time(),'\n', sep=''))

#Extract WAIC up to here
WAIC<-Cgsg.MCMC$getWAIC()

saveRDS(WAIC, file=paste0('./result_chain1/waic_run',i,'.rds', sep=''))
cat(paste0('saved WAIC 5/5, run ',i,' at ',Sys.time(),'\n', sep=''))

rm(stateList, chain1, chain1_logProbs, WAIC)
end.chain_r1=Sys.time()



cat(paste0('finished run ',i,' at ',Sys.time(),' in ',(end.chain_r1-start.chain), '\n', sep=''))

}

# start.chain=Sys.time()
# #Run 2.1
# Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# end.chain=Sys.time()
# 
# #Extract and save samples
# chain1_run2<-as.matrix(Cgsg.MCMC$mvSamples)
# saveRDS(chain1_run2,file = './result_chain1/chain1_run2.rds')
# 
# chain1_run2_age<-as.matrix(Cgsg.MCMC$mvSamples2)
# saveRDS(chain1_run2_age,file = './result_chain1/chain1_run2_age.rds')
# 
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run2.rds')
# 
# rm(stateList, WAIC)
# 
# #Run 2.2
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run2<-rbind(chain1_run2,samples)
# saveRDS(chain1_run2,file = './result_chain1/chain1_run2.rds')
# 
# 
# chain1_run2_age<-rbind(chain1_run2_age,samples2)
# saveRDS(chain1_run2_age,file = './result_chain1/chain1_run2_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run2.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# 
# #Run 2.3
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run2<-rbind(chain1_run2,samples)
# saveRDS(chain1_run2,file = './result_chain1/chain1_run2.rds')
# 
# 
# chain1_run2_age<-rbind(chain1_run2_age,samples2)
# saveRDS(chain1_run2_age,file = './result_chain1/chain1_run2_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run2.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 2.4
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run2<-rbind(chain1_run2,samples)
# saveRDS(chain1_run2,file = './result_chain1/chain1_run2.rds')
# 
# 
# chain1_run2_age<-rbind(chain1_run2_age,samples2)
# saveRDS(chain1_run2_age,file = './result_chain1/chain1_run2_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run2.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 2.5
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run2<-rbind(chain1_run2,samples)
# saveRDS(chain1_run2,file = './result_chain1/chain1_run2.rds')
# 
# 
# chain1_run2_age<-rbind(chain1_run2_age,samples2)
# saveRDS(chain1_run2_age,file = './result_chain1/chain1_run2_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run2.rds')
# 
# rm(samples, samples2, stateList, chain1_run2, chain1_run2_age, WAIC)
# end.chain_r1=Sys.time()
# (time.chain<-end.chain_r1-start.chain)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# start.chain=Sys.time()
# #Run 3.1
# Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# end.chain=Sys.time()
# 
# #Extract and save samples
# chain1_run3<-as.matrix(Cgsg.MCMC$mvSamples)
# saveRDS(chain1_run3,file = './result_chain1/chain1_run3.rds')
# 
# chain1_run3_age<-as.matrix(Cgsg.MCMC$mvSamples2)
# saveRDS(chain1_run3_age,file = './result_chain1/chain1_run3_age.rds')
# 
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run3.rds')
# 
# rm(stateList, WAIC)
# 
# #Run 3.2
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run3<-rbind(chain1_run3,samples)
# saveRDS(chain1_run3,file = './result_chain1/chain1_run3.rds')
# 
# 
# chain1_run3_age<-rbind(chain1_run3_age,samples2)
# saveRDS(chain1_run3_age,file = './result_chain1/chain1_run3_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run3.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# 
# #Run 3.3
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run3<-rbind(chain1_run3,samples)
# saveRDS(chain1_run3,file = './result_chain1/chain1_run3.rds')
# 
# 
# chain1_run3_age<-rbind(chain1_run3_age,samples2)
# saveRDS(chain1_run3_age,file = './result_chain1/chain1_run3_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run3.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 3.4
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run3<-rbind(chain1_run3,samples)
# saveRDS(chain1_run3,file = './result_chain1/chain1_run3.rds')
# 
# 
# chain1_run3_age<-rbind(chain1_run3_age,samples2)
# saveRDS(chain1_run3_age,file = './result_chain1/chain1_run3_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run3.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 3.5
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run3<-rbind(chain1_run3,samples)
# saveRDS(chain1_run3,file = './result_chain1/chain1_run3.rds')
# 
# 
# chain1_run3_age<-rbind(chain1_run3_age,samples2)
# saveRDS(chain1_run3_age,file = './result_chain1/chain1_run3_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run3.rds')
# 
# rm(samples, samples2, stateList, chain1_run3, chain1_run3_age, WAIC)
# end.chain_r1=Sys.time()
# (time.chain<-end.chain_r1-start.chain)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# start.chain=Sys.time()
# #Run 4.1
# Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# end.chain=Sys.time()
# 
# #Extract and save samples
# chain1_run4<-as.matrix(Cgsg.MCMC$mvSamples)
# saveRDS(chain1_run4,file = './result_chain1/chain1_run4.rds')
# 
# chain1_run4_age<-as.matrix(Cgsg.MCMC$mvSamples2)
# saveRDS(chain1_run4_age,file = './result_chain1/chain1_run4_age.rds')
# 
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run4.rds')
# 
# rm(stateList, WAIC)
# 
# #Run 4.2
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run4<-rbind(chain1_run4,samples)
# saveRDS(chain1_run4,file = './result_chain1/chain1_run4.rds')
# 
# 
# chain1_run4_age<-rbind(chain1_run4_age,samples2)
# saveRDS(chain1_run4_age,file = './result_chain1/chain1_run4_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run4.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# 
# #Run 4.3
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run4<-rbind(chain1_run4,samples)
# saveRDS(chain1_run4,file = './result_chain1/chain1_run4.rds')
# 
# 
# chain1_run4_age<-rbind(chain1_run4_age,samples2)
# saveRDS(chain1_run4_age,file = './result_chain1/chain1_run4_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run4.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# 
# #Run 4.4
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run4<-rbind(chain1_run4,samples)
# saveRDS(chain1_run4,file = './result_chain1/chain1_run4.rds')
# 
# 
# chain1_run4_age<-rbind(chain1_run4_age,samples2)
# saveRDS(chain1_run4_age,file = './result_chain1/chain1_run4_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run4.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 4.5
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run4<-rbind(chain1_run4,samples)
# saveRDS(chain1_run4,file = './result_chain1/chain1_run4.rds')
# 
# 
# chain1_run4_age<-rbind(chain1_run4_age,samples2)
# saveRDS(chain1_run4_age,file = './result_chain1/chain1_run4_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run4.rds')
# 
# rm(samples, samples2, stateList, chain1_run4, chain1_run4_age, WAIC)
# end.chain_r1=Sys.time()
# (time.chain<-end.chain_r1-start.chain)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# start.chain=Sys.time()
# #Run 5.1
# Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# end.chain=Sys.time()
# 
# #Extract and save samples
# chain1_run5<-as.matrix(Cgsg.MCMC$mvSamples)
# saveRDS(chain1_run5,file = './result_chain1/chain1_run5.rds')
# 
# chain1_run5_age<-as.matrix(Cgsg.MCMC$mvSamples2)
# saveRDS(chain1_run5_age,file = './result_chain1/chain1_run5_age.rds')
# 
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run5.rds')
# 
# rm(stateList, WAIC)
# 
# #Run 5.2
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run5<-rbind(chain1_run5,samples)
# saveRDS(chain1_run5,file = './result_chain1/chain1_run5.rds')
# 
# 
# chain1_run5_age<-rbind(chain1_run5_age,samples2)
# saveRDS(chain1_run5_age,file = './result_chain1/chain1_run5_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run5.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 5.3
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run5<-rbind(chain1_run5,samples)
# saveRDS(chain1_run5,file = './result_chain1/chain1_run5.rds')
# 
# 
# chain1_run5_age<-rbind(chain1_run5_age,samples2)
# saveRDS(chain1_run5_age,file = './result_chain1/chain1_run5_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run5.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 5.4
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run5<-rbind(chain1_run5,samples)
# saveRDS(chain1_run5,file = './result_chain1/chain1_run5.rds')
# 
# 
# chain1_run5_age<-rbind(chain1_run5_age,samples2)
# saveRDS(chain1_run5_age,file = './result_chain1/chain1_run5_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run5.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 5.5
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run5<-rbind(chain1_run5,samples)
# saveRDS(chain1_run5,file = './result_chain1/chain1_run5.rds')
# 
# 
# chain1_run5_age<-rbind(chain1_run5_age,samples2)
# saveRDS(chain1_run5_age,file = './result_chain1/chain1_run5_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run5.rds')
# 
# rm(samples, samples2, stateList, chain1_run5, chain1_run5_age, WAIC)
# end.chain_r1=Sys.time()
# (time.chain<-end.chain_r1-start.chain)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# start.chain=Sys.time()
# #Run 6.1
# Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# end.chain=Sys.time()
# 
# #Extract and save samples
# chain1_run6<-as.matrix(Cgsg.MCMC$mvSamples)
# saveRDS(chain1_run6,file = './result_chain1/chain1_run6.rds')
# 
# chain1_run6_age<-as.matrix(Cgsg.MCMC$mvSamples2)
# saveRDS(chain1_run6_age,file = './result_chain1/chain1_run6_age.rds')
# 
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run6.rds')
# 
# rm(stateList, WAIC)
# 
# #Run 6.2
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run6<-rbind(chain1_run6,samples)
# saveRDS(chain1_run6,file = './result_chain1/chain1_run6.rds')
# 
# 
# chain1_run6_age<-rbind(chain1_run6_age,samples2)
# saveRDS(chain1_run6_age,file = './result_chain1/chain1_run6_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run6.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# 
# #Run 6.3
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run6<-rbind(chain1_run6,samples)
# saveRDS(chain1_run6,file = './result_chain1/chain1_run6.rds')
# 
# 
# chain1_run6_age<-rbind(chain1_run6_age,samples2)
# saveRDS(chain1_run6_age,file = './result_chain1/chain1_run6_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run6.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 6.4
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run6<-rbind(chain1_run6,samples)
# saveRDS(chain1_run6,file = './result_chain1/chain1_run6.rds')
# 
# 
# chain1_run6_age<-rbind(chain1_run6_age,samples2)
# saveRDS(chain1_run6_age,file = './result_chain1/chain1_run6_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run6.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 6.5
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run6<-rbind(chain1_run6,samples)
# saveRDS(chain1_run6,file = './result_chain1/chain1_run6.rds')
# 
# 
# chain1_run6_age<-rbind(chain1_run6_age,samples2)
# saveRDS(chain1_run6_age,file = './result_chain1/chain1_run6_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run6.rds')
# 
# rm(samples, samples2, stateList, chain1_run6, chain1_run6_age, WAIC)
# end.chain_r1=Sys.time()
# (time.chain<-end.chain_r1-start.chain)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# start.chain=Sys.time()
# #Run 7.1
# Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# end.chain=Sys.time()
# 
# #Extract and save samples
# chain1_run7<-as.matrix(Cgsg.MCMC$mvSamples)
# saveRDS(chain1_run7,file = './result_chain1/chain1_run7.rds')
# 
# chain1_run7_age<-as.matrix(Cgsg.MCMC$mvSamples2)
# saveRDS(chain1_run7_age,file = './result_chain1/chain1_run7_age.rds')
# 
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run7.rds')
# 
# rm(samples, WAIC)
# 
# #Run 7.2
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run7<-rbind(chain1_run7,samples)
# saveRDS(chain1_run7,file = './result_chain1/chain1_run7.rds')
# 
# 
# chain1_run7_age<-rbind(chain1_run7_age,samples2)
# saveRDS(chain1_run7_age,file = './result_chain1/chain1_run7_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run7.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# 
# #Run 7.3
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run7<-rbind(chain1_run7,samples)
# saveRDS(chain1_run7,file = './result_chain1/chain1_run7.rds')
# 
# 
# chain1_run7_age<-rbind(chain1_run7_age,samples2)
# saveRDS(chain1_run7_age,file = './result_chain1/chain1_run7_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run7.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 7.4
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run7<-rbind(chain1_run7,samples)
# saveRDS(chain1_run7,file = './result_chain1/chain1_run7.rds')
# 
# 
# chain1_run7_age<-rbind(chain1_run7_age,samples2)
# saveRDS(chain1_run7_age,file = './result_chain1/chain1_run7_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run7.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 7.5
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run7<-rbind(chain1_run7,samples)
# saveRDS(chain1_run7,file = './result_chain1/chain1_run7.rds')
# 
# 
# chain1_run7_age<-rbind(chain1_run7_age,samples2)
# saveRDS(chain1_run7_age,file = './result_chain1/chain1_run7_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run7.rds')
# 
# rm(samples, samples2, stateList, chain1_run7, chain1_run7_age, WAIC)
# end.chain_r1=Sys.time()
# (time.chain<-end.chain_r1-start.chain)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# start.chain=Sys.time()
# #Run 8.1
# Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# end.chain=Sys.time()
# 
# #Extract and save samples
# chain1_run8<-as.matrix(Cgsg.MCMC$mvSamples)
# saveRDS(chain1_run8,file = './result_chain1/chain1_run8.rds')
# 
# chain1_run8_age<-as.matrix(Cgsg.MCMC$mvSamples2)
# saveRDS(chain1_run8_age,file = './result_chain1/chain1_run8_age.rds')
# 
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run8.rds')
# 
# rm(samples, WAIC)
# 
# 
# #Run 8.2
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run8<-rbind(chain1_run8,samples)
# saveRDS(chain1_run8,file = './result_chain1/chain1_run8.rds')
# 
# 
# chain1_run8_age<-rbind(chain1_run8_age,samples2)
# saveRDS(chain1_run8_age,file = './result_chain1/chain1_run8_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run8.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# 
# #Run 8.3
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run8<-rbind(chain1_run8,samples)
# saveRDS(chain1_run8,file = './result_chain1/chain1_run8.rds')
# 
# 
# chain1_run8_age<-rbind(chain1_run8_age,samples2)
# saveRDS(chain1_run8_age,file = './result_chain1/chain1_run8_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run8.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 8.4
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run8<-rbind(chain1_run8,samples)
# saveRDS(chain1_run8,file = './result_chain1/chain1_run8.rds')
# 
# 
# chain1_run8_age<-rbind(chain1_run8_age,samples2)
# saveRDS(chain1_run8_age,file = './result_chain1/chain1_run8_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run8.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 8.5
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run8<-rbind(chain1_run8,samples)
# saveRDS(chain1_run8,file = './result_chain1/chain1_run8.rds')
# 
# 
# chain1_run8_age<-rbind(chain1_run8_age,samples2)
# saveRDS(chain1_run8_age,file = './result_chain1/chain1_run8_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run8.rds')
# 
# rm(samples, samples2, stateList, chain1_run8, chain1_run8_age, WAIC)
# end.chain_r1=Sys.time()
# (time.chain<-end.chain_r1-start.chain)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# start.chain=Sys.time()
# #Run 9.1
# Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# end.chain=Sys.time()
# 
# #Extract and save samples
# chain1_run9<-as.matrix(Cgsg.MCMC$mvSamples)
# saveRDS(chain1_run9,file = './result_chain1/chain1_run9.rds')
# 
# chain1_run9_age<-as.matrix(Cgsg.MCMC$mvSamples2)
# saveRDS(chain1_run9_age,file = './result_chain1/chain1_run9_age.rds')
# 
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run9.rds')
# 
# rm(samples, WAIC)
# 
# #Run 9.2
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run9<-rbind(chain1_run9,samples)
# saveRDS(chain1_run9,file = './result_chain1/chain1_run9.rds')
# 
# 
# chain1_run9_age<-rbind(chain1_run9_age,samples2)
# saveRDS(chain1_run9_age,file = './result_chain1/chain1_run9_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run9.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# 
# #Run 9.3
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run9<-rbind(chain1_run9,samples)
# saveRDS(chain1_run9,file = './result_chain1/chain1_run9.rds')
# 
# 
# chain1_run9_age<-rbind(chain1_run9_age,samples2)
# saveRDS(chain1_run9_age,file = './result_chain1/chain1_run9_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run9.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 9.4
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run9<-rbind(chain1_run9,samples)
# saveRDS(chain1_run9,file = './result_chain1/chain1_run9.rds')
# 
# 
# chain1_run9_age<-rbind(chain1_run9_age,samples2)
# saveRDS(chain1_run9_age,file = './result_chain1/chain1_run9_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run9.rds')
# 
# rm(samples, samples2, stateList, WAIC)
# 
# #Run 9.5
# Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
# 
# #Extract and save samples
# samples<-as.matrix(Cgsg.MCMC$mvSamples)
# samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
# 
# #Extract and save samples
# chain1_run9<-rbind(chain1_run9,samples)
# saveRDS(chain1_run9,file = './result_chain1/chain1_run9.rds')
# 
# 
# chain1_run9_age<-rbind(chain1_run9_age,samples2)
# saveRDS(chain1_run9_age,file = './result_chain1/chain1_run9_age.rds')
# 
# #Extract model and MCMC state and save them
# stateList <- list(modelState = getModelState(C.gsg.mod),
#                   mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
#                   rs = .Random.seed)
# 
# saveRDS(stateList, file = './result_chain1/model_state.rds')
# 
# 
# #Extract WAIC up to here
# WAIC<-Cgsg.MCMC$getWAIC()
# 
# saveRDS(WAIC, file='./result_chain1/waic_run9.rds')
# 
# rm(samples, samples2, stateList, chain1_run9, chain1_run9_age, WAIC)
# end.chain_r1=Sys.time()
# (time.chain<-end.chain_r1-start.chain)
