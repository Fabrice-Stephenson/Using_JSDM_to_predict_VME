##===========================  S3 FIT MODEL OBJECTS       =====================================##

##------------------------------  AUTHORS: Fabrice Stephenson based on code from Ovaskainen and Abrego 2020
##------------------------------  Project : Fisheries New Zealand funded project: ZBD201901
##------------------------------  Start date : 01/03/2019
##------------------------------  End date : 05/07/2023

##=============================================================================================================##
library(Hmsc)

# THIS SCRIPT FITS HMSC MODELS BASED ON CODE FROM THE THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# Set working directory 
wd <- ""
setwd(wd)
load(file = "models/unfitted models") #models, modelnames

# Number of samples to iteratively run
samples_list <- c(5,200,200)
# thinning parameter increasing until convergence
thin_list <- c(1,1,10)
# Number of chains (which can be run in parralell) 
nChains <- 5

# models[[1]]$XFormula

for(Lst in 1:length(samples_list)){
  thin <- thin_list[Lst]
  samples <- samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  nm <- length(models)
  for (model in 1:nm) {
    start <- Sys.time()
    print(paste0("model = ",modelnames[model]))
    m <- models[[model]]
    m <- sampleMcmc(m, samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains, nParallel = nChains) #nParallel
    models[[model]] <- m
    end <- Sys.time()
    print(end - start)
  }
  filename <- paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")
  save(models,modelnames,file=filename)
}
