##===========================  S8 COMPUTE SPATIAL PREDICTIONS               =====================================##

##------------------------------  AUTHORS: Fabrice Stephenson based on code from Ovaskainen and Abrego 2020
##------------------------------  Project : Fisheries New Zealand funded project: ZBD201901
##------------------------------  Start date : 01/03/2019
##------------------------------  End date : 05/07/2023

##=============================================================================================================##
library(Hmsc); library(ggplot2); library(raster)

# THIS SCRIPT COMPUTES SPATIAL PREDICTIONS FOR FITTED HMSC MODELS BASED ON CODE FROM THE THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# Set working directory 
wd <- ""
setwd(wd)

localDir <- "."

# Only calculating for the "final" model
nChains <- 5
samples <- 200
thin <- 10
filename = paste("models_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),
                 ".Rdata",sep = "")
load(filename)
nm = length(models)

# Load dataframe which represents environemtnal variables used as predictors in the JSDM 
# for each 1 km grid of the study area
load("ZBDpred_1km.Rdata")
crs <- CRS("+proj=merc +lat_ts=-41 +lon_0=100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

ZBDpred_1km <- na.omit(ZBDpred_1km)

# Because the model predictions are memory intensive, the area is split into smaller spaces
# for prediction, averaged, and then reassembled. 

# SUBSET TO STUDY AREA
test_area <- ZBDpred_1km[ZBDpred_1km$y <(-3600000) & ZBDpred_1km$y > (-5700000),]
# r <- rasterFromXYZ(data.frame(x = test_area[,1],
#                               y =test_area[,2],
#                               z =test_area[,3]),
#                    crs = crs)
# plot(r)

# xy <- as.matrix(cbind(test_area$x, test_area$y))
xy <- as.matrix(cbind(test_area$x, test_area$y))
colnames(xy)=c("X","Y")

save(xy, file = "xy.Rdata")

imp.vars <- colnames(models[[1]]$XData)

test_area <- (test_area[,imp.vars])

for (i in 1:ncol(test_area)){
  r <- rasterFromXYZ(data.frame(x =xy[,1],
                                y =xy[,2],
                                z =test_area[,i]),
                     crs = crs)
  names(r) <- colnames(test_area[i])
  if (i ==1){Pred.stack <- r}else{Pred.stack <- stack(Pred.stack,r)}
}

Pred.stack
plot(Pred.stack[[1]])

# subset the data into 30 * 30 areas (10800 chunks)
XdataN <- na.omit(as.data.frame(Pred.stack, xy = T))
min(XdataN$x)
max(XdataN$x)
min(XdataN$y)
max(XdataN$y)

x.split <- seq(min(XdataN$x),max(XdataN$x),length = 30)
y.split <- seq(min(XdataN$y),max(XdataN$y),length = 30)

counter <- 0
squares <- list()
for (i in 1:(length(x.split)-1)){
  for (j in 1:(length(y.split)-1)){
  counter <- counter +1
  squares[[counter]] <- data.frame(as.matrix(XdataN[XdataN$x >= x.split[i] & XdataN$x < x.split[i+1]
                   & XdataN$y >= y.split[j] & XdataN$y < y.split[j+1],]))
  }
}

# remove all dataframes with 0 observations (i.e. NAs)
squares <- squares[sapply(squares, nrow)>0]

# check the number of cells on average for calc of time
mean.length <- c()
for (i in 1:length(squares)){mean.length[i] <- (nrow(squares[[i]]))}
sum(mean.length)/length(squares)

# SETUP PARRALLELISATION 
library(parallel);library(foreach);library(bigstatsr)
# detectCores()

packages <- c("Hmsc")

multiResultClass <- function(mean.pred=NULL, UC.pred=NULL) 
{me <- list(mean.pred = mean.pred, UC.pred = UC.pred)
## Set the name for the class
class(me) <- append(class(me),"multiResultClass")
return(me)
}

cl <- parallel::makeCluster(10) # number of cores to use - reduce this to avoid overburdening the server
doParallel::registerDoParallel(cl)

# length(squares)
start <- Sys.time()
# PARALLEL LOOPS
tmp3 <- foreach(i = 1:length(squares)) %dopar% {
  lapply(packages, require, character.only = TRUE) # load other packages
  
  xy <- as.matrix(cbind(squares[[i]]$x, squares[[i]]$y))
  colnames(xy)=c("X","Y")
  XdataN <- squares[[i]][,-c(1:2)]
  
  sp.grad <- prepareGradient(models[[1]], XDataNew = XdataN,
                             sDataNew = list(Transect_ID = xy))
  predY <- predict(models[[1]], Gradient=sp.grad, expected = F ,predictETAMean = T)
  
  a.PA <- simplify2array(predY) # turn to arrray for averaging
  EpredY <- apply(a.PA, c(1,2), mean) # mean value 
  EpredY.UC <- apply(a.PA, c(1,2), var) # uncert
  
  result <- multiResultClass()
  
  result$mean.pred <- cbind(xy,EpredY)
  result$UC.pred <- cbind(xy,EpredY.UC)
  
  return(result)
}

parallel::stopCluster(cl)
# insert serial backend, otherwise error in repetetive tasks
registerDoSEQ()

end <- Sys.time()
end - start


for(j in 1:length(tmp3)){
  if (j == 1){predY_PA <- tmp3[[j]]$mean.pred
  predY_PA.UC <- tmp3[[j]]$UC.pred
  } else {
    predY_PA <- rbind(predY_PA, tmp3[[j]]$mean.pred)
    predY_PA.UC <- rbind(predY_PA.UC, tmp3[[j]]$UC.pred)
  }
}

save(predY_PA,predY_PA.UC, file = "predY_PA.Rdata")

# Abundance 
cl <- parallel::makeCluster(10) # number of cores to use - reduce this to avoid overburdening the server
doParallel::registerDoParallel(cl)

# length(squares)
start <- Sys.time()
# PARALLEL LOOPS
tmp3 <- foreach(i = 1:length(squares)) %dopar% {
  lapply(packages, require, character.only = TRUE) # load other packages
  
  xy <- as.matrix(cbind(squares[[i]]$x, squares[[i]]$y))
  colnames(xy)=c("X","Y")
  XdataN <- squares[[i]][,-c(1:2)]
  
  sp.grad <- prepareGradient(models[[2]], XDataNew = XdataN,
                             sDataNew = list(Transect_ID = xy))
  predY <- predict(models[[2]], Gradient=sp.grad, expected = F ,predictETAMean = T)
  
  a.PA <- simplify2array(predY) # turn to arrray for averaging
  EpredY <- apply(a.PA, c(1,2), mean) # mean value 
  EpredY.UC <- apply(a.PA, c(1,2), var) # uncert
  
  result <- multiResultClass()
  
  result$mean.pred <- cbind(xy,EpredY)
  result$UC.pred <- cbind(xy,EpredY.UC)
  
  return(result)
}

parallel::stopCluster(cl)
# insert serial backend, otherwise error in repetetive tasks
registerDoSEQ()

end <- Sys.time()
end - start


for(j in 1:length(tmp3)){
  if (j == 1){PredY_Abun <- tmp3[[j]]$mean.pred
  PredY_Abun.UC <- tmp3[[j]]$UC.pred
  } else {
    PredY_Abun <- rbind(PredY_Abun, tmp3[[j]]$mean.pred)
    PredY_Abun.UC <- rbind(PredY_Abun.UC, tmp3[[j]]$UC.pred)
  }
}

save(PredY_Abun, PredY_Abun.UC, file = "predY_Abun.Rdata")


####--------------    DISTRIBUTION MAP            -------------------------------#### 
library(raster)
crs <- CRS("+proj=merc +lat_ts=-41 +lon_0=100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

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

wd <- ""
setwd(wd)
# Mean PA + uncert
load("figures/predY_PA.Rdata")
load("figures/predY_Abun.Rdata")
load("figures/xy.Rdata")

mask <- raster("Mask_1500.tif")
plot(mask)

for (i in 1:length(sp.code)){
  Mpred.r <- rasterFromXYZ(data.frame(x =predY_PA[,1],
                                      y =predY_PA[,2],
                                      z =predY_PA[,sp.code[i]]),
                           crs = crs)
  SDpred.r <- rasterFromXYZ(data.frame(x =predY_PA.UC[,1],
                                       y =predY_PA.UC[,2],
                                       z =predY_PA.UC[,sp.code[i]]),
                            crs = crs)
  
  # abundance
  Mpred2.r <- rasterFromXYZ(data.frame(x =PredY_Abun[,1],
                                       y =PredY_Abun[,2],
                                       z =PredY_Abun[,sp.code[i]]),
                            crs = crs)
  SDpred2.r <- rasterFromXYZ(data.frame(x =PredY_Abun.UC[,1],
                                        y =PredY_Abun.UC[,2],
                                        z =PredY_Abun.UC[,sp.code[i]]),
                             crs = crs)
  M.Dens.r <- Mpred.r * Mpred2.r
  SDDens.r <- SDpred.r * SDpred2.r
  
  writeRaster(Mpred.r, filename = paste0("rasters/",sp.code[i],"_PA.tif"))
  writeRaster(SDpred.r, filename = paste0("rasters/",sp.code[i],"_PA_SD.tif"))
  writeRaster(M.Dens.r, filename = paste0("rasters/",sp.code[i],"_DENS.tif"))
  writeRaster(SDDens.r, filename = paste0("rasters/",sp.code[i],"_DENS_SD.tif"))
}

# SPECIES RICHNESS
getS <- function(x){return(rowSums(x))}
# aS <- lapply(predY_PA,getS)
# head(predY_PA)

S <- as.data.frame(cbind(predY_PA[,1:2],rowSums(predY_PA[,3:ncol(predY_PA)])))
S.UC <- as.data.frame(cbind(predY_PA.UC[,1:2],rowSums(predY_PA.UC[,3:ncol(predY_PA.UC)])))

Richness <- rasterFromXYZ(data.frame(x =S[,1],
                                     y =S[,2],
                                     z =S[,3]),
                          crs = crs)
Richness.UC <- rasterFromXYZ(data.frame(x =S.UC[,1],
                                        y =S.UC[,2],
                                        z =S.UC[,3]),
                             crs = crs)
plot(Richness.UC)
writeRaster(Richness, filename = paste0("rasters/","Richness.tif"))
writeRaster(Richness.SD, filename = paste0("rasters/","Richness.SD.tif"))
