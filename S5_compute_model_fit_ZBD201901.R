##===========================  S5 COMPUTE MODEL FIT                =====================================##

##------------------------------  AUTHORS: Fabrice Stephenson based on code from Ovaskainen and Abrego 2020
##------------------------------  Project : Fisheries New Zealand funded project: ZBD201901
##------------------------------  Start date : 01/03/2019
##------------------------------  End date : 05/07/2023

##=============================================================================================================##
library(Hmsc)

# THIS SCRIPT COMPUTES MODEL FIT METRICS FOR FITTED HMSC MODELS BASED ON CODE FROM THE THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# Set working directory 
wd <- ""
setwd(wd)

# Only calculating for the "final" model
nChains <- 5
thin <- 10
samples <- 200
print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
filename_in <- paste("models/models_thin_", as.character(thin),
                    "_samples_", as.character(samples),
                    "_chains_",as.character(nChains),
                    ".Rdata",sep = "")

load(file = filename_in) #models, modelnames
nm <- length(models)

MF <- list()
MFCV <- list()
WAIC <- list()

# Loop through both Presence/absence model and abundance model
# Performance metrics on model fit with data, and on 10-fold cross-validation
for(model in 1:nm){
  print(paste0("model = ",as.character(model)))
  start <- Sys.time()
  m <- models[[model]]
  preds <- computePredictedValues(m)
  MF[[model]] <- evaluateModelFit(hM=m, predY=preds)
  partition <- createPartition(m, nfolds = 10)
  preds <- computePredictedValues(m,partition=partition, nParallel = 5) # make parrallel
                                 # , partition.sp = 1:m$ns, mcmcStep = 10) #sp partition
  MFCV[[model]] <- evaluateModelFit(hM=m, predY=preds)
  WAIC[[model]] <- computeWAIC(m)
  end <- Sys.time()
  print(end - start)
}

filename_out <- paste("models/MF_thin_conditional_", as.character(thin),
                                            "_samples_", as.character(samples),
                                            "_chains_",as.character(nChains),
                                            ".Rdata",sep = "")
save(MF,MFCV,WAIC,modelnames,file = filename_out)
