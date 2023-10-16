##===========================  S2 DEFINE MODEL OBJECTS       =================================================##

##------------------------------  AUTHORS: Fabrice Stephenson based on code from Ovaskainen and Abrego 2020
##------------------------------  Project : Fisheries New Zealand funded project: ZBD201901
##------------------------------  Start date : 01/03/2019
##------------------------------  End date : 05/07/2023

##=============================================================================================================##
library(Hmsc)

# THIS SCRIPT CONSTRUCTS HMSC MODELS BASED ON CODE FROM THE THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# create the subdirectories "data" and "models" to the main directory

# Set working directory 
wd <- ""
setwd(wd)

localDir <- "."
data.directory <- file.path(localDir, "data")
model.directory <- file.path(localDir, "models")
load(file = file.path(data.directory, "allData.R")) # S, X, Y & Tr Dataframes

##### ---------------  EXPLORATION OF SPECIES DATA     ---------------------#####

# check for absent (0) or ubitquitous species (1)
range(colMeans(Y>0))
min(colMeans(Y>0))

# Create a dataframe of rare taxa
# IMPORTANT - we could split out REEF into rare taxa here
raretaxa <- which(colSums(Y>0) < 15)
length(raretaxa)

# remove rare taxa from the full species data - in this case nothing removed
# Y <- Y[,-raretaxa]

hist(colMeans(Y > 0 ), main = "Prevalence", breaks = 20)
hist(log(Y[Y > 0]), main = "Log abundance conditional on presence",breaks = 20)
# this looks more or less normal 
hist(rowSums(Y>0), main = "Richness across samples", breaks = 20) # species richness across samples
# also looks fairly normal -interesting 

# Transect ID that correponds to the identity number of the transect.
head(S)
S$Transect_ID <- droplevels(S$Transect_ID)
# The data contains also the x- and y-coordinates of the routes. We store these as the xy-matrix to be able to fit a spatial model
xy <- data.frame(cbind(S$X,S$Y))

rownames(xy)<- S$Transect_ID
colnames(xy) <- c("X","Y")
par(mfrow=c(1,1))
plot(xy, asp=1) # show the map (NB., equal aspect ratio in the map)

# Let us then look at the environmental and spatial data
head(X)
names(X)

# the following also show the types (factor, num) of variables:
str(X)
# and a full summary
summary(X)
plot(X)

#####-----------         DEFINING THE HMSC MODEL                           -----------------#######
# the model includes both the environmental covariates as well as the spatial random effect of the transect ID 

studyDesign <- data.frame(Transect_ID = as.factor(S$Transect_ID))

# Defining the random level object. 
rL <- HmscRandomLevel(sData=xy)

# Note that the row names of xy correspond to the units of the studyDesign. This is necessary to make
# Hmsc understand how the units of the random effect (rows of xy) correspond to the sampling units
# (rows of studyDesign, Y and XData). 

# creation of the HMSC formula
XFormula = ~  Bathy + footprint_NZR + Salinity_NZR + 
  Tempres_NZR + prof_curv_NZR + TidalCurr_NZR + epp_me_NZR
     
# FINAL Y DATA FOR BOTH PA AND ABUNDANCE MODELS
# given the 0s and the abundance data we will need to make a hurdle model 
# to account for P/A and separately for abundance conditional on presence

Ypa <- 1*(Y>0)
Yabu <- Y
Yabu[Y == 0] <- NA
Yabu <- log(Yabu)
# hist(Yabu)

# DEFINING MODEL
# Presence/absence model
m1 <- Hmsc(Y=Ypa, XData = X, XFormula=XFormula,
         # phyloTree = phyloTree, 
         # TrData = Tr, TrFormula = TrFormula,
         distr="probit", studyDesign=studyDesign,
         ranLevels=list(Transect_ID=rL))

# Abundance (conditional on presence) model
m2 <-  Hmsc(Y=Yabu, YScale = T, XData = X, XFormula=XFormula,
            # phyloTree = phyloTree, 
            # TrData = Tr, TrFormula = TrFormula,
            distr="normal", studyDesign=studyDesign,
            ranLevels=list(Transect_ID=rL))

models <- list(m1, m2)
modelnames <- c("Presence_absence", "Abundance_COP") 
save(models, modelnames, file = file.path(model.directory, "unfitted models"))