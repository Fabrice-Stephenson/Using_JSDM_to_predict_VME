##===========================  S9 MAPPING SPATIAL PREDICTIONS & VME PREDICTIONS         =====================================##

##------------------------------  AUTHORS: Fabrice Stephenson based on code from Ovaskainen and Abrego 2020
##------------------------------  Project : Extension of FNZ project - funded by 
##------------------------------            Sustainable Seas National Science Challenge Project 3.2
##------------------------------  Start date : 01/03/2019
##------------------------------  End date : 05/07/2023

##=============================================================================================================##
library(Hmsc); library(ggplot2); library(raster); library(rgdal);library(sp);library(spdep)

# THIS SCRIPT MAPS SPATIAL PREDICTIONS & ESTIMATES VME PREDICTIONS as described in 
# Stephenson et al. (in review) 'Using joint species distribution modelling to predict distributions of seafloor 
# taxa and identify vulnerable marine ecosystems in New Zealand waters'

# Set working directory 
wd <- ""
setwd(wd)

# MAPPING TAXA              -----------------------------------------------------------------------------------

#### PROJECTIONS
merc41proj <- CRS("+init=epsg:3994")
wgs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
MPIproj <- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

#### LOAD JSDM PREDICTIONS
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

# extensions
ext <- c("_PA.tif", "_PA_SD.tif", "_Dens.tif", "_Dens_SD.tif")

PA <- c()
PA_SD <- c()
Dens <- c()
Dens_SD <- c()

#PA
for (i in 1:length(sp.code)){
  setwd(paste(wd, "/rasters/", sep = ""))
  r <- raster(paste(sp.code[[i]], ext[[1]], sep = ""))
  names(r) <- sp.code[[i]]
  if (i == 1){
    PA <- stack(r)
  } else {
    PA <- stack(PA, r)
  }
}

#PA SD
for (i in 1:length(sp.code)){
  setwd(paste(wd, "/rasters/", sep = ""))
  r <- raster(paste(sp.code[[i]], ext[[2]], sep = ""))
  names(r) <- sp.code[[i]]
  if (i == 1){
    PA_SD <- stack(r)
  } else {
    PA_SD <- stack(PA_SD, r)
  }
}

#DENS
for (i in 1:length(sp.code)){
  setwd(paste(wd, "/rasters/", sep = ""))
  r <- raster(paste(sp.code[[i]], ext[[3]], sep = ""))
  r <- exp(r)
  r <- r * 1000
  names(r) <- sp.code[[i]]
  if (i == 1){
    Dens <- stack(r)
  } else {
    Dens <- stack(Dens, r)
  }
}
#DENS SD
for (i in 1:length(sp.code)){
  setwd(paste(wd, "/rasters/", sep = ""))
  r <- raster(paste(sp.code[[i]], ext[[4]], sep = ""))
  r <- r * 1000
  # r <- exp(r)
  names(r) <- sp.code[[i]]
  if (i == 1){
    Dens_SD <- stack(r)
  } else {
    Dens_SD <- stack(Dens_SD, r)
  }
}

# plot(Dens[["Goniocorella.dumosa_Dens"]])
plot(Dens[[16]])

#### BASEMAP LAYERS
setwd("")
nzcoast <- readOGR(dsn = "NZ_coastline_Albers.shp")
nzcoast <- spTransform(nzcoast, merc41proj)
nzcoast.r <- raster::rasterize(nzcoast, PA[[1]], "rast")
nzcoast <- nzcoast.r

# plot(nzcoast, add = T, col="black", border = NA)

#### LABELS FOR MAPS
xlb<-spTransform(SpatialPoints(cbind(c(160,170,180,190),c(-42,-42,-42,-42)),proj4string=CRS("+proj=longlat")), merc41proj)
ylb<-spTransform(SpatialPoints(cbind(c(176,176,176,176,176),c(-40,-45,-50,-55,-60)),proj4string=CRS("+proj=longlat")), merc41proj)
legend_pos <- spTransform(SpatialPoints(cbind(c(190),c(-56)),proj4string=CRS("+proj=longlat")), merc41proj)
# legend_pos@coords
Fig_ext <- spTransform(SpatialPoints(cbind(c(160,185),c(-42,-57)),proj4string=CRS("+proj=longlat")), merc41proj)

Lx<-c("160","170","180","170")
Ly<-c("40","45","50","55","60")
Lx2<-parse(text = paste(Lx, "*degree", sep = ""))
Ly2<-parse(text = paste(Ly, "*degree", sep = ""))

#### PLOTTING AND SAVING MAPS TO PDF
setwd(paste(wd, "/Figures"))
pdf(file = "PA.DensTaxa.pdf")
for(i in 1:nlayers(PA)) {
  if (any(i == c(1:38,40:42,44:62, 64:67))){
    #Just the updated models
    # png(filename = paste(sp.code[i],".png", sep=""), width = 20, 
    #     height = 20, units = 'cm', bg = "white", res = 1000)
    par(mar=c(2,2,1.5,0),mgp=c(1.8,0.8,0))
    par(mfrow=c(2,2))
    
    # PA
    brk.PA <- seq(0,1,0.1)
    col.PA <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
    r <- PA[[i]]
    r <- extend(r, c(100,100))
    plot(r, colNA=NA, col=col.PA, breaks=brk.PA, legend=FALSE,
         bty="n", axes = F, box = T, asp = 1,
         xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
         main =paste0(sp.code[i], " - Prob. of occ.")) 
    plot(nzcoast,add=TRUE,col="black", border = NA, legend=FALSE)
    axis(side=1, at=xlb$coords.x1,labels=Lx2, cex.axis = 1.0)
    axis(side=2, pos = xlb$coords.x1[1], at=ylb$coords.x2,labels=Ly2, cex.axis = 1.0)
    plot(r,legend.only=TRUE,col=col.PA, breaks=brk.PA,
         alpha=1,legend.width=1,legend.shrink=0.75,
         smallplot=c(0.62,0.67,0.13,0.43), main = "test"
         # legend.args=list(text='Prob Occ', font=2, line=2.5, cex=0.8)
    )
    
    # PA.SD
    brk.PA.SD <- round(seq(0,0.5, length.out = 9),2)
    col.PA.SD <- RColorBrewer::brewer.pal(9,"Reds")
    r <- PA_SD[[i]]
    r <- extend(r, c(100,100))
    plot(r, colNA=NA, col=col.PA.SD, breaks=brk.PA.SD, legend=FALSE,
         bty="n", axes = F, box = T, asp = 1,
         xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
         main =paste0("SD"))
    plot(nzcoast,add=TRUE,col="black", border = NA, legend=FALSE)
    axis(side=1, at=xlb$coords.x1,labels=Lx2, cex.axis = 1.0)
    axis(side=2, pos = xlb$coords.x1[1], at=ylb$coords.x2,labels=Ly2, cex.axis = 1.0)
    plot(r,legend.only=TRUE,col=col.PA.SD, breaks=brk.PA.SD,
         alpha=1,legend.width=1,legend.shrink=0.75,
         smallplot=c(0.62,0.67,0.13,0.43), main = "test"
         # legend.args=list(text='Prob Occ', font=2, line=2.5, cex=0.8)
    )
    
    # Dens
    q <- quantile(Dens[[i]], probs = c(0.98), na.rm = T)
    brk.Dens <- c(seq(min(na.omit(values(Dens[[i]]))), q,length.out = 9),max(na.omit(values(Dens[[i]]))))
    col.Dens <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
    r <- Dens[[i]]
    r <- extend(r, c(100,100))
    plot(r, colNA=NA, col=col.Dens, breaks=brk.Dens, legend=FALSE,
         bty="n", axes = F, box = T, asp = 1,
         xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
         main =paste0(sp.code[i], " - Density")) 
    plot(nzcoast,add=TRUE,col="black", border = NA, legend=FALSE)
    axis(side=1, at=xlb$coords.x1,labels=Lx2, cex.axis = 1.0)
    axis(side=2, pos = xlb$coords.x1[1], at=ylb$coords.x2,labels=Ly2, cex.axis = 1.0)
    
    r.leg <- clamp(r, 0, q)
    brk.Dens <- unique(round(c(seq(min(na.omit(values(Dens[[i]]))), q,length.out = 10)),0))
    col.Dens <- rev(RColorBrewer::brewer.pal(length(brk.Dens),"RdBu"))
    plot(r.leg,legend.only=TRUE,col=col.Dens, breaks=brk.Dens,2,
         alpha=1,legend.width=1,legend.shrink=1.0,
         smallplot=c(0.62,0.67,0.13,0.43)
         # legend.args=list(text='Prob Occ', font=2, line=2.5, cex=0.8)
    )
    
    # Dens.SD
    q <- quantile(Dens_SD[[i]], probs = c(0.98), na.rm = T)
    brk.Dens.SD <- round(c(seq(min(na.omit(values(Dens_SD[[i]]))), q,length.out = 8),max(na.omit(values(Dens_SD[[i]])))),0)
    col.Dens.SD <- RColorBrewer::brewer.pal(9,"Reds")
    r <- Dens_SD[[i]]
    r <- extend(r, c(100,100))
    plot(r, colNA=NA, col=col.Dens.SD, breaks=brk.Dens.SD, legend=FALSE,
         bty="n", axes = F, box = T, asp = 1,
         xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
         main =paste0("SD")) 
    plot(nzcoast,add=TRUE,col="black", border = NA, legend=FALSE)
    axis(side=1, at=xlb$coords.x1,labels=Lx2, cex.axis = 1.0)
    axis(side=2, pos = xlb$coords.x1[1], at=ylb$coords.x2,labels=Ly2, cex.axis = 1.0)
    
    r.leg <- clamp(r, 0, q)
    brk.Dens.SD <- unique(round(c(seq(min(na.omit(values(Dens_SD[[i]]))), q, length.out = 9)),0))
    col.Dens.SD <- rev(RColorBrewer::brewer.pal(length(brk.Dens.SD),"Reds"))
    plot(r.leg,legend.only=TRUE,col=col.Dens.SD, breaks=brk.Dens.SD,2,
         alpha=1,legend.width=1,legend.shrink=1.0,
         smallplot=c(0.62,0.67,0.13,0.43))
    print(paste("Finished iteration: ", i, "of 66." ))
  } else {
    par(mar=c(2,2,1.5,0),mgp=c(1.8,0.8,0))
    par(mfrow=c(2,2))
    
    # PA
    brk.PA <- seq(0,1,0.1)
    col.PA <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
    r <- PA[[i]]
    r <- extend(r, c(100,100))
    plot(r, colNA=NA, col=col.PA, breaks=brk.PA, legend=FALSE,
         bty="n", axes = F, box = T, asp = 1,
         xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
         main =paste0(sp.code[i], " - Prob. of occ.")) 
    plot(nzcoast,add=TRUE,col="black", border = NA, legend=FALSE)
    axis(side=1, at=xlb$coords.x1,labels=Lx2, cex.axis = 1.0)
    axis(side=2, pos = xlb$coords.x1[1], at=ylb$coords.x2,labels=Ly2, cex.axis = 1.0)
    plot(r,legend.only=TRUE,col=col.PA, breaks=brk.PA,
         alpha=1,legend.width=1,legend.shrink=0.75,
         smallplot=c(0.62,0.67,0.13,0.43), main = "test"
         # legend.args=list(text='Prob Occ', font=2, line=2.5, cex=0.8)
    )
    
    # PA.SD
    brk.PA.SD <- round(seq(0,0.5, length.out = 9),2)
    col.PA.SD <- RColorBrewer::brewer.pal(9,"Reds")
    r <- PA_SD[[i]]
    r <- extend(r, c(100,100))
    plot(r, colNA=NA, col=col.PA.SD, breaks=brk.PA.SD, legend=FALSE,
         bty="n", axes = F, box = T, asp = 1,
         xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
         main =paste0("SD"))
    plot(nzcoast,add=TRUE,col="black", border = NA, legend=FALSE)
    axis(side=1, at=xlb$coords.x1,labels=Lx2, cex.axis = 1.0)
    axis(side=2, pos = xlb$coords.x1[1], at=ylb$coords.x2,labels=Ly2, cex.axis = 1.0)
    plot(r,legend.only=TRUE,col=col.PA.SD, breaks=brk.PA.SD,
         alpha=1,legend.width=1,legend.shrink=0.75,
         smallplot=c(0.62,0.67,0.13,0.43), main = "test"
         # legend.args=list(text='Prob Occ', font=2, line=2.5, cex=0.8)
    )
    
    # Dens
    q <- quantile(Dens[[i]], probs = c(0.98), na.rm = T)
    brk.Dens <- c(seq(min(na.omit(values(Dens[[i]]))), q,length.out = 9))
    col.Dens <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
    r <- Dens[[i]]
    r <- extend(r, c(100,100))
    plot(r, colNA=NA, col=col.Dens, breaks=brk.Dens, legend=FALSE,
         bty="n", axes = F, box = T, asp = 1,
         xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
         main =paste0(sp.code[i], " - Density")) 
    plot(nzcoast,add=TRUE,col="black", border = NA, legend=FALSE)
    axis(side=1, at=xlb$coords.x1,labels=Lx2, cex.axis = 1.0)
    axis(side=2, pos = xlb$coords.x1[1], at=ylb$coords.x2,labels=Ly2, cex.axis = 1.0)
    
    r.leg <- clamp(r, 0, q)
    brk.Dens <- unique(round(c(seq(min(na.omit(values(Dens[[i]]))), q,length.out = 10)),0))
    col.Dens <- rev(RColorBrewer::brewer.pal(length(brk.Dens),"RdBu"))
    plot(r.leg,legend.only=TRUE,col=col.Dens, breaks=brk.Dens,2,
         alpha=1,legend.width=1,legend.shrink=1.0,
         smallplot=c(0.62,0.67,0.13,0.43)
         # legend.args=list(text='Prob Occ', font=2, line=2.5, cex=0.8)
    )
    
    # Dens.SD
    q <- quantile(Dens_SD[[i]], probs = c(0.98), na.rm = T)
    brk.Dens.SD <- round(c(seq(min(na.omit(values(Dens_SD[[i]]))), q,length.out = 8)),0)
    col.Dens.SD <- RColorBrewer::brewer.pal(9,"Reds")
    r <- Dens_SD[[i]]
    r <- extend(r, c(100,100))
    plot(r, colNA=NA, col=col.Dens.SD, breaks=brk.Dens.SD, legend=FALSE,
         bty="n", axes = F, box = T, asp = 1,
         xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
         main =paste0("SD")) 
    plot(nzcoast,add=TRUE,col="black", border = NA, legend=FALSE)
    axis(side=1, at=xlb$coords.x1,labels=Lx2, cex.axis = 1.0)
    axis(side=2, pos = xlb$coords.x1[1], at=ylb$coords.x2,labels=Ly2, cex.axis = 1.0)
    
    r.leg <- clamp(r, 0, q)
    brk.Dens.SD <- unique(round(c(seq(min(na.omit(values(Dens_SD[[i]]))), q, length.out = 9)),0))
    col.Dens.SD <- rev(RColorBrewer::brewer.pal(length(brk.Dens.SD),"Reds"))
    plot(r.leg,legend.only=TRUE,col=col.Dens.SD, breaks=brk.Dens.SD,2,
         alpha=1,legend.width=1,legend.shrink=1.0,
         smallplot=c(0.62,0.67,0.13,0.43))
    print(paste("Finished iteration: ", i, "of 66." ))
  }
}
dev.off()

# MAPPING VME               -----------------------------------------------------------------------------------

#### TAXA, MORPHOTYPE, VULNERABILTY SCORE
# load info on vulnerability for each VME taxa (not the others)
setwd("")
VME.score <- read.csv('VME_vulnerabilty.csv', sep = ";")

VME.taxa <- VME.score$VME
# VME.taxa %in% sp.code # checking that all VME taxa have same names as saved in models

#### ABUNDANCE WEIGHTED VULNERABILTY 

for (i in 1:length(VME.taxa)){
  r <- Dens[[VME.taxa[i]]] # extract spatial prediction based on VME taxa name
  r <- r * VME.score[VME.score$VME == VME.taxa[i],2] # multiply spatial pred by vulnerabilty score based on taxa name
  if (i == 1){
    Abu.Vul <- r
  } else {
    Abu.Vul <- Abu.Vul + r
  }
}
# dev.new()

# core areas
Abu.Vul.top5 <- Abu.Vul
Abu.Vul.top5[values(Abu.Vul.top5) < quantile(Abu.Vul, c(0.95))] <- 0
Abu.Vul.top5[values(Abu.Vul.top5) >= quantile(Abu.Vul, c(0.95))] <- 1

plot(Abu.Vul)
plot(Abu.Vul.top5)

#### RICHNESS WEIGHTED VULNERABILTY / RICHNESS
for (i in 1:length(VME.taxa)){
  r <- PA[[VME.taxa[i]]] # extract spatial prediction based on VME taxa name
  r <- r * VME.score[VME.score$VME == VME.taxa[i],2] # multiply spatial pred by vulnerabilty score based on taxa name
  if (i == 1){
    Rich.Vul <- r
  } else {
    Rich.Vul <- Rich.Vul + r
  }
}


Rich.Vul.top5 <- Rich.Vul
Rich.Vul.top5[values(Rich.Vul.top5) < quantile(Rich.Vul, c(0.95))] <- 0
Rich.Vul.top5[values(Rich.Vul.top5) >= quantile(Rich.Vul, c(0.95))] <- 2

plot(Rich.Vul)
plot(Rich.Vul.top5)

# OVERLAP OF METRICS
VME.cmb.5 <- Rich.Vul.top5 + Abu.Vul.top5

#### UNCERTAINTY WEIGHTED VULNERABILTY 
for (i in 1:length(VME.taxa)){
  r <- Dens_SD[[VME.taxa[i]]] # extract spatial prediction based on VME taxa name
  # r <- r * VME.score[VME.score$VME == VME.taxa[i],2] # multiply spatial pred by vulnerabilty score based on taxa name
  if (i == 1){
    UC.Vul <- r
  } else {
    UC.Vul <- UC.Vul + r
  }
}
dev.new()
plot(UC.Vul)

CV.Vul <- UC.Vul/Abu.Vul
plot(CV.Vul)

UC.Vul.top10 <- CV.Vul
UC.Vul.top10[values(UC.Vul.top10) < 0.33] <- 0

# dev.new()
plot(UC.Vul.top10)

# PLOT METRICS IN HIGH RES PDF TO ALLOW ZOOM IN

# 1. Abundance weighted vulnerability for each VME taxa (30)
setwd(paste(wd, "/Figures"))
pdf(file = "Abu.Vul.pdf")
par(mar=c(2,2,1.5,0),mgp=c(1.8,0.8,0))

for(i in 1:length(VME.taxa)) {
    par(mar=c(2,2,1.5,0),mgp=c(1.8,0.8,0))
    
    r <- Dens[[VME.taxa[i]]] # extract spatial prediction based on VME taxa name
    r <- r * VME.score[VME.score$VME == VME.taxa[i],2]
    # plot(r)
    
    # Dens
    q <- quantile(r, probs = c(0.98), na.rm = T)
    brk.Dens <- c(seq(min(na.omit(values(r))), q,length.out = 9),max(na.omit(values(r))))
    col.Dens <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
    r <- extend(r, c(100,100))
    plot(r, colNA=NA, col=col.Dens, breaks=brk.Dens, legend=FALSE,
         bty="n", axes = F, box = T, asp = 1,
         xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
         main =paste0(VME.taxa[i], " - Vulnerabilty Index")) 
    plot(nzcoast,add=TRUE,col="black", border = NA, legend = F)
    axis(side=1, at=xlb$coords.x1,labels=Lx2, cex.axis = 1.0)
    axis(side=2, pos = xlb$coords.x1[1], at=ylb$coords.x2,labels=Ly2, cex.axis = 1.0)
    
    r.leg <- clamp(r, 0, q)
    brk.Dens <- unique(round(c(seq(min(na.omit(values(r))), q,length.out = 10)),2))
    col.Dens <- rev(RColorBrewer::brewer.pal(length(brk.Dens),"RdBu"))
    plot(r.leg,legend.only=TRUE,col=col.Dens, breaks=brk.Dens,2,
         alpha=1,legend.width=1,legend.shrink=1.0,
         smallplot=c(0.62,0.67,0.13,0.43))
    print(paste("Finished iteration: ", i, "of 30." ))
}
dev.off()
  
# 2.	Abundance weighted VME index 
setwd("")
pdf(file = "Abu.Ind.pdf")
par(mar=c(2,2,1.5,0),mgp=c(1.8,0.8,0))
par(mfrow=c(2,2))

plot(Abu.Vul, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="Abundance weighted VME Ind")
plot(nzcoast,add=TRUE,col="black", border = NA)
plot(Abu.Vul.top10, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="VME Ind (Abu) - Top 10%")
plot(nzcoast,add=TRUE,col="black", border = NA)
plot(Abu.Vul.top5, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="VME Ind (Abu) - Top 5%")
plot(nzcoast,add=TRUE,col="black", border = NA)
plot(Abu.Vul.Jenk, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="VME Ind (Abu) - Jenks")
plot(nzcoast,add=TRUE,col="black", border = NA)
dev.off()

# 3.	Richness weighted VME index 
setwd("")
pdf(file = "Rich.Ind.pdf")
par(mar=c(2,2,1.5,0),mgp=c(1.8,0.8,0))
par(mfrow=c(2,2))
plot(Rich.Vul, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="Richness weighted VME Ind")
plot(nzcoast,add=TRUE,col="black", border = NA)
plot(Rich.Vul.top10, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="VME Ind (Rich) - Top 10%")
plot(nzcoast,add=TRUE,col="black", border = NA)
plot(Rich.Vul.top5, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="VME Ind (Rich) - Top 5%")
plot(nzcoast,add=TRUE,col="black", border = NA)
plot(Rich.Vul.Jenk, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="VME Ind (Rich) - Jenks")
plot(nzcoast,add=TRUE,col="black", border = NA)
dev.off()

# 4.	Intersection of VME indices 
setwd("")
pdf(file = "VME.Ind.pdf")
par(mar=c(2,2,1.5,0),mgp=c(1.8,0.8,0))
par(mfrow=c(2,2))

brk<- c(0,0.5,1,2,3)
col <- c("gray92", "darkorange", "darkorchid1", "deepskyblue")

plot(VME.cmb.10, col=col, breaks=brk, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="Overlap 10% (1: abu; 2 = rich; 3 = both)")
plot(nzcoast,add=TRUE,col="black", border = NA)
plot(VME.cmb.5, col=col, breaks=brk, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="Overlap 5 (1: abu; 2 = rich; 3 = both)")
plot(nzcoast,add=TRUE,col="black", border = NA)
plot(VME.cmb.Jenk, col=col, breaks=brk, colNA=NA, legend=T,
     bty="n", axes = F, box = T, asp = 1,
     xlim = c(5048111, 7300000), ylim = c(-5800000, -3750000),
     main ="Overlap Jenk (1: abu; 2 = rich; 3 = both)")
plot(nzcoast,add=TRUE,col="black", border = NA)

dev.off()

#### SAVE RASTERS
setwd("")
writeRaster(Abu.Vul, file = "Abu.Ind.tif", overwrite = T)
writeRaster(Rich.Vul, file = "Rich.Ind.tif", overwrite = T)
writeRaster(Abu.Vul.top5, file = "Abu.Ind.Top5.tif", overwrite = T)
writeRaster(Rich.Vul.top5, file = "Rich.Ind.Top5.tif", overwrite = T)
writeRaster(VME.cmb.5, file = "VME.Ind.Top5.tif", overwrite = T)
writeRaster(UC.Vul, file = "VME.Ind.UC.tif", overwrite = T)
writeRaster(UC.Vul.top10, file = "VME.Ind.UC.Thresh.tif", overwrite = T)


