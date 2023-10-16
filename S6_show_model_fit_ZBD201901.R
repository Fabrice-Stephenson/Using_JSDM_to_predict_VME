##===========================  S6 VISUALISE MODEL FITS              =====================================##

##------------------------------  AUTHORS: Fabrice Stephenson based on code from Ovaskainen and Abrego 2020
##------------------------------  Project : Fisheries New Zealand funded project: ZBD201901
##------------------------------  Start date : 01/03/2019
##------------------------------  End date : 05/07/2023

##=============================================================================================================##
library(Hmsc); library(writexl)

# THIS SCRIPT VISUALISES MODEL FITS FOR FITTED HMSC MODELS BASED ON CODE FROM THE THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# Set working directory 
wd <- ""
setwd(wd)

# Only calculating for the "final" model
thin <- 10
samples <- 200
nChains <- 5

filename <- paste("models/MF_thin_conditional_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),
                 ".Rdata",sep = "")

load(filename)
nm <- length(MF)
filename <- paste("figures/model_fit.pdf")
pdf(file = filename)

sp.code <- c("Isididae", "Bivalvia", "Anemones", "Serolidae", "Buccinidae", "Alcyonacea", 
             "Gastropoda", "Polychaeta", "Hexactinellida", "Anthomastus.sp.", "Pennatulacea", 
             "Antipatharia", "Ascidiacea", "Holothuroidea", "Asteroidea", "Demospongiae", 
             "Gorgonocephalidae", "Brachyura", "Barnacles", "Ophiuroidea", "Octopoda", 
             "Bryozoa", "Brachiopoda", "Brisingidae", "Euechinoidea", "Primnoidae", 
             "Stylasteridae", "Caridea", "Ceriantharia", "Galatheidae.Chirostylidae", 
             "Gorgonacea", "Cidaroidea", "Cladorhizidae", "Pycnogonida", "Corallimorpharia", 
             "Scleractinia", "Crinoidea..motile.", "Crinoidea..stalked.", "Crustacean..lobster.",
             "Paguridae", "Metanephrops.challengeri", "Caryophylliidae", "Stephanocyathus.sp.", 
             "Dermechinus.horridus", "Nudibranchia", "Echinothurioida", "Echiura", 
             "Enypniastes.eximia", "Epizoanthidae", "Flabellum", "Xenophyophoroidea", "Ranellidae", 
             "Goniocorella.dumosa", "Hyalascus.n..sp.", "Hydrozoa", "Kophobelemnon.sp.", 
             "Worm.indet.", "Spatangoida", "Psolidae", "Hyalinoecia.sp.", "Radicipes.sp.", 
             "Scaphopoda", "Siphonophore", "Taiaroa.tauhou", "Telesto.sp.", "Volutidae", 
             "Zoanthidea")

excel.list <- list()

# Display as graphs and as excel file 

for(j in 1:nm){
  cMF <- MF[[j]]
  cMFCV <- MFCV[[j]]
  if(!is.null(cMF$TjurR2)){
    plot(cMF$TjurR2,cMFCV$TjurR2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": Tjur R2"))
    abline(0,1)
    abline(v=0)
    abline(h=0)
    
    # create excel file for reference
    me <- data.frame(cbind(cMF$TjurR2,cMFCV$TjurR2))
    me <- data.frame(cbind(sp.code,me))
    colnames(me) <- c("Taxa","TjurR2: Explanatory ", "TjurR2: Predictive")
    excel.list <- c(excel.list, "TjurR2" = list(me))
  }
  if(!is.null(cMF$R2)){
    plot(cMF$R2,cMFCV$R2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(modelnames[[j]],", thin = ",as.character(thin),", 
                     samples = ",as.character(samples),": R2"))
    abline(0,1)
    abline(v=0)
    abline(h=0)
    
    # create excel file for reference
    me <- data.frame(cbind(cMF$R2,cMFCV$R2))
    me <- data.frame(cbind(sp.code,me))
    colnames(me) <- c("Taxa","R2: Explanatory ", "R2: Predictive")
    rownames(me) <- sp.code
    excel.list <- c(excel.list, "R2" = list(me))
  }
  if(!is.null(cMF$AUC)){
    plot(cMF$AUC,cMFCV$AUC,xlim=c(0,1),ylim=c(0,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": AUC"))
    abline(0,1)
    abline(v=0.5)
    abline(h=0.5)
    
    # create excel file for reference
    me <- data.frame(cbind(cMF$AUC,cMFCV$AUC))
    me <- data.frame(cbind(sp.code,me))
    colnames(me) <- c("Taxa","AUC: Explanatory ", "AUC: Predictive")
    rownames(me) <- sp.code
    excel.list <- c(excel.list, "AUC" = list(me))
  }
  
  filename <- paste0("figures/cross_val_",modelnames[j],".xlsx")
  writexl::write_xlsx(excel.list,path = filename)
}
dev.off()
