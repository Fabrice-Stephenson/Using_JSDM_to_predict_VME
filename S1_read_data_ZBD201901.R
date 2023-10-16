##===========================  S1 READ DATA AND PREP MODEL OBJECTS       =====================================##

##------------------------------  AUTHORS: Fabrice Stephenson based on code from Ovaskainen and Abrego 2020
##------------------------------  Project : Fisheries New Zealand funded project: ZBD201901
##------------------------------  Start date : 01/03/2019
##------------------------------  End date : 05/07/2023

##=============================================================================================================##

##=====================             RAW DATA              ======================================================
# Provide an SXY file.
#set working directory
wd <- ""
setwd(wd)
load("biol_AGG_clean.Rdata")
colnames(DF_biol)
DF <- DF_biol; rm(DF_biol) 

# Data frame with columns x, Y coordinates (Mercator -41 projection), Transect ID, taxa and environmental predictors
# for each site sampled (rows)
# head(DF$Transect_ID)

# only keep those with more than 10 samples
sp <- DF[,8:81]
sp.code <- colnames(sp[,which(colSums(sp>0) > 9)]); rm(sp)

# All variables available
# "Bathy_NZR"     "BPI_broad_NZR" "cbpm_me_NZR"   "Dissox_NZR"    "DOM_NZR"      
# "Dynoc_NZR"     "epp_me_NZR"    "epp_min_NZR"   "Flux_new"      "footprint_NZR"
# "poc_NZR"       "prof_curv_NZR" "Salinity_NZR"  "sil_NZR"       "sstgrad_NZR"  
# "std15_NZR"     "Tempres_NZR"   "TidalCurr_NZR"
# imp.var <- colnames(DF[,28:ncol(DF)])
imp.var <- c("Bathy", "Flux_new", "footprint_NZR",
             "Salinity_NZR","Tempres_NZR", "prof_curv_NZR","TidalCurr_NZR","epp_me_NZR")

####--------       CO-LINEARITY OF ENVIRONMENTAL PREDICTORS      ----------------------####
cordf <- DF[,imp.var] # columns with predictors for the whole dataset 
cut <- 4 # VIF cutoff, 3 - 4 is usually the arbitrarily selected cutoff 

repeat {
  corvif<-usdm::vif(cordf)  # Where cordf is a dataframe of variables
  x<-max(corvif[,2])      # You can play around with different values of x from Yesson et al.,2015
  if (x<cut){break}
  else {cordf<-cordf[,-which.max(corvif[,2])]}
}

corvif

# CORRELATION OF PREDICTOR VARIABLES
cordf <- na.omit(cordf)
cordf<-cor(cordf)
head(round(cordf,2))

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(cordf)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# save figure
setwd(paste(wd, "/figures", sep ="")) # file directory for saving files
devEMF::emf(file = "Corr.plot.emf", emfPlus = FALSE)
corrplot::corrplot(cordf, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)
# rerun all the graph code to "print" it to the blank EMF
dev.off()

# final variables 
imp.var <- colnames(cordf)

####--------       SETTING UP THE FILE STRUCTURE              ----------------------####
# The files TP and P are optional, so indicate with TRUE/FALSE if they are included or not
is.TP = F
is.P = F

# READING IN SXY: study design (S) and/or covariates (X) and species data (Y) 
SXY <- DF; rm(DF)
SXY <- na.omit(SXY)
S <- SXY[,2:4]
X <- SXY[,imp.var]
Y <- SXY[,sp.code]

# Check that the data looks as it should
View(S)
View(X)
View(Y)

# check that community data are numeric and have finite numbers. 
if (is.numeric(as.matrix(Y)) || is.logical(as.matrix(Y)) && is.finite(sum(Y, na.rm=TRUE))) {
    print("Y looks OK")
} else {
	print("Y should be numeric and have finite values")	}
# Check that the stydy design data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(S))) {
  print("S has NA values - not allowed for")
} else {
  print("S looks ok")	}
# Check that the covariate data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(X))) {
  print("X has NA values - not allowed for")
} else {
  print("X looks ok")	}

setwd(wd)
# READING IN TP: traits (T) and/or phylogenetic information in table format (P)
if(is.TP){
  # Read in the species names as rownames, not as a column of the matrix
  TP <- read.csv("TP.csv", stringsAsFactors=TRUE,row.names = 1)
  # The script below checks if the species names in TP are identical and in the same order as in Y
  # If the script prints "species names in TP and SXY match", you are ok.
  # If it says that they do not match, you need to modify the files so that they match 
  if(all(rownames(TP)==colnames(Y))) {
    print("species names in TP and SXY match")
  } else{
    print("species names in TP and SXY do not match")
  }
  # Modify the next two lines to split your TP file to components that relate to
  # Tr: species traits (note that T is a reserved word in R and that's why we use Tr)
  # P: phylogenetic information given by taxonomical levels, e.g. order, family, genus, species
  # If you don't have trait data, indicate this by Tr=NULL. 
  # If TP does not have phylogenetic data (because you don't have such data at all, or because
  # it is given in tree-format, like is the case in this example), indicate this with P=NULL 
  Tr = TP[,1:ncol(TP)]
  P = NULL
  # Check that the data looks as it should!
  View(Tr)
  View(P)
  # Check that the Tr data do not have missing values (they are allowed for Y but not S, X, P or Tr)
  if (any(is.na(Tr))) {
    print("Tr has NA values - not allowed for")
  } else {
    print("Tr looks ok")	}
  # Check that the phylogenetic/taxonomic data do not have missing values (they are allowed for Y but not S, X, P or Tr)
  if (any(is.na(P))) {
    print("P has NA values - not allowed for")
  } else {
    print("P looks ok")	}
}

# READING IN P: phylogenetic information in tree format (P)
# we use ape package for trees, and P.tre must be in a format that ape understands
if(is.P){
  # Read in the phylogenetic tree using read.tree from ape
  library(ape)
  P = read.tree("P.tre")
  # When you look at P (e.g. write P and press enter),
  # you should see that it is a phylogenetic tree which
  # is rooted and includes branch lengths and tip labels
  # The script below checks if the species names in P are identical (but not necessarily in the same order) as in Y
  # If the script prints "species names in P and SXY match", you are ok.
  # If it says that they do not match, you need to modify the files so that they match 
  if(all(sort(P$tip.label) == sort(colnames(Y)))){
    print("species names in P and SXY match")
  } else{
    print("species names in P and SXY do not match")
  }
  # Check that the data looks as it should!
  plot(P, cex=0.5)
}

####--------       FINAL FILES SAVED FOR NEXT SECTION              ----------------------####
setwd(paste(wd, "/data", sep = ""))
save(S, X, Y, file = "allData.R")
