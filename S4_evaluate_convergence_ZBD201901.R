##===========================  S4 EVALUATE MODEL CONVERGENCE       =====================================##

##------------------------------  AUTHORS: Fabrice Stephenson based on code from Ovaskainen and Abrego 2020
##------------------------------  Project : Fisheries New Zealand funded project: ZBD201901
##------------------------------  Start date : 01/03/2019
##------------------------------  End date : 05/07/2023

##=============================================================================================================##
library(Hmsc)

# THIS SCRIPT EVALUATES MODEL CONVERGENCE OF FITTED HMSC MODELS BASED ON CODE FROM THE THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# Set working directory 
wd <- ""
setwd(wd)

#include in samples_list and thin_list only those models that you have actually fitted!
samples_list <- c(5,200,200)
thin_list <- c(1,1,10)
nst <- length(thin_list)
nChains <- 5

ma <- NULL
na <- NULL
for (Lst in 1:nst) {
  thin <- thin_list[Lst]
  samples <- samples_list[Lst]
  
  
  filename <- paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),".Rdata",sep = "")
  load(filename)
  nm <- length(models)
  for(j in 1:nm){
    mpost <- convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta <- gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    tmp <- summary(psrf.beta)
    if(is.null(ma)){
      ma <- psrf.beta[,1]
      na <- paste0(as.character(thin),",",as.character(samples))
    } else {
      ma <- cbind(ma,psrf.beta[,1])
      if(j==1){
        na <- c(na,paste0(as.character(thin),",",as.character(samples)))
      } else {
        na <- c(na,"")
      }
    }
  }
}

pdf(file=paste("figures/MCMC_convergence.pdf"))
par(mfrow=c(2,1))
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
dev.off()

