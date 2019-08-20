# Maps -------------------------------------------------------------------
#Modern sites
library(maps)
library(mapdata)
library(raster)
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Plots/Maps")
modern_pollen<- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv", row.names=1)
m<-map('world', col = 'gray90', fill = T, border = 'gray60',
    bg = 'white',xlim=c(min(modern_pollen$Long),max(modern_pollen$Long)),
    ylim=c(min(modern_pollen$Lat),max(modern_pollen$Lat)),mar = c(5,2, par("mar")[3], 0.1))

xat <- pretty(m$range[1:2])
xlab <- parse(text=degreeLabelsEW(xat))
yat <- pretty(m$range[3:4])
ylab <- parse(text=degreeLabelsNS(yat))
box()
axis(1, at=xat, labels=xlab)
axis(2, las=TRUE, at=yat, labels=ylab)
axis(3, at=xat, labels=xlab)
axis(4, las=TRUE, at=yat, labels=ylab)

points(modern_pollen$Long,modern_pollen$Lat,col="gray30",pch=16,cex=0.5)

#Fossil sites
h1 <- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/h1.csv")
h1<-h1[,-1]
h1<-h1[order(h1$longitude),]
m<-map('world', col = 'gray90', fill = T, border = 'gray60',
    bg = 'white', xlim=c(-10,4),ylim=c(35,44),mar = c(0, 0, par("mar")[3], 0.1))

xat <- pretty(m$range[1:2])
xlab <- parse(text=degreeLabelsEW(xat))
yat <- pretty(m$range[3:4])
ylab <- parse(text=degreeLabelsNS(yat))
box()
axis(1, at=xat, labels=xlab)
axis(2, las=TRUE, at=yat, labels=ylab)
axis(3, at=xat, labels=xlab)
axis(4, las=TRUE, at=yat, labels=ylab)

points(h1$longitude,h1$latitude,col="black",pch=16,cex=1)

text(h1$longitude+0.15,h1$latitude,h1$name,cex= 0.75, pos=2)
text(h1$longitude+0.3,h1$latitude-0.1,h1$altitude,cex= 0.7, pos=3)

# Core -------------------------------------------------------------
rm(list=ls())
install.packages("rioja")
install.packages("ggplot2")
install.packages("maps")

modern_pollen<- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv", row.names=1)
taxaColMin<-which(colnames(modern_pollen)=="Abies")
taxaColMax<-which(colnames(modern_pollen)=="Zygophyllaceae")
Holocene <- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Holocene.csv", row.names=1)
core<-Holocene[,-c(1:3)]

source('C:/Users/ml4418/Desktop/Master Project/Script/Tolerance weighted WA-PLS/Define functions.R', encoding = 'UTF-8')

# WA-PLS ------------------------------------------------------------------
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Input data")
component<-3
pred_core<-WAPLS_core(modern_pollen,taxaColMin,taxaColMax,component,core)
pred_core<-cbind.data.frame(Holocene[,c(1:3)],pred_core)

write.csv(pred_core,"pred_core.csv")

# Tolerance weighted WA-PLS -----------------------------------------------
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Input data")
component<-2
P<-read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Optimum and Tolerance of each taxa.csv", row.names=1)
pred_t_core<-Tolerance_weighted_WAPLS_core(modern_pollen,taxaColMin,taxaColMax,P,component,core)
pred_t_core<-cbind.data.frame(Holocene[,c(1:3)],pred_t_core)

write.csv(pred_t_core,"pred_t_core.csv")

# Standard error by DKL ---------------------------------------------------
SE_gdd_DKL<-Tolerance_weighted_WAPLS_SE_DKL(core,P$t1)
SE_alpha_DKL<-Tolerance_weighted_WAPLS_SE_DKL(core,P$t2)
SE_Tmin_DKL<-Tolerance_weighted_WAPLS_SE_DKL(core,P$t3)

SE_DKL<-cbind.data.frame(SE_gdd_DKL,SE_alpha_DKL,SE_Tmin_DKL)
pred_t_core<-cbind.data.frame(pred_t_core,SE_DKL)

write.csv(SE_DKL,"SE_DKL.csv")
write.csv(pred_t_core,"pred_t_core.csv")

# Standard error by sampling the training set 1000 times----------------------------------------------------------
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Input data")

# Tolerance weighted WAPLS Standard Error by sampling the trainin --------
component<-2
train<- modern_pollen[,taxaColMin:taxaColMax]
core<-Holocene[,-c(1:3)]

replace1000<-Tolerance_weighted_WAPLS_SE_replace1000(train, modern_pollen$gdd, t(P$sd1), t(P$t1),n=1000, component, core = core) 
y2_gdd<- apply(replace1000, 1, sd, na.rm=TRUE)                  # Get SDs
recm_gdd<- apply(replace1000, 1, mean, na.rm=TRUE)   # Get means

replace1000<-Tolerance_weighted_WAPLS_SE_replace1000(train, modern_pollen$alpha, t(P$sd2), t(P$t2),n=1000, component, core) 
y2_alpha<- apply(replace1000, 1, sd, na.rm=TRUE)                  # Get SDs
recm_alpha<- apply(replace1000, 1, mean, na.rm=TRUE)   # Get means

replace1000<-Tolerance_weighted_WAPLS_SE_replace1000(train,modern_pollen$Tmin, t(P$sd3), t(P$t3),n=1000, component, core) 
y2_Tmin<- apply(replace1000, 1, sd, na.rm=TRUE)                  # Get SDs
recm_Tmin<- apply(replace1000, 1, mean, na.rm=TRUE)   # Get means

SE_replace1000<-cbind.data.frame(y2_gdd,y2_alpha,y2_Tmin,recm_gdd,recm_alpha,recm_Tmin)
names(SE_replace1000)<-c("SE_gdd_replace1000","SE_alpha_replace1000","SE_Tmin_replace1000","gdd_replace1000","alpha_replace1000","Tmin_replace1000")
write.csv(SE_replace1000,"SE_replace1000.csv")

pred_t_core<-cbind.data.frame(pred_t_core,SE_replace1000[,c("SE_gdd_replace1000","SE_alpha_replace1000","SE_Tmin_replace1000")])
write.csv(pred_t_core,"pred_t_core.csv")

# Insolation --------------------------------------------------------------
insol_site<-data.frame(matrix(nrow=nrow(pred_t_core),ncol=4))
for(j in 1:nrow(pred_t_core)){
  insol_site[j,1:4]<-insol_summer_winter[which(insol_summer_winter$Site==as.character(pred_t_core[j,"Site"])&insol_summer_winter$Age==as.numeric(pred_t_core[j,"Age.cal.BP"])),]
}
pred_t_core$insol_summer<-insol_site[,3]
pred_t_core$insol_winter<-insol_site[,4]
write.csv(pred_t_core,"pred_t_core.csv")

#######################################################################################################

# Plot the climates changes during the Holocene ---------------------------
#Get the modern climate of the fossil sites
library(readxl)
library(readr)
h1 <- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/h1.csv")
h1<-h1[,-1]
h1<-h1[order(h1$longitude),]

#Output for thesis
climate_fossil<-h1[,c("name","gdd","pred_t_gdd","alpha","pred_t_alpha","Tmin","pred_t_Tmin")]
climate_fossil[,c("gdd","pred_t_gdd")]<-round(climate_fossil[,c("gdd","pred_t_gdd")],digits=0)
climate_fossil[,c("alpha","pred_t_alpha")]<-round(climate_fossil[,c("alpha","pred_t_alpha")],digits=3)
climate_fossil[,c("Tmin","pred_t_Tmin")]<-round(climate_fossil[,c("Tmin","pred_t_Tmin")],digits=2)
write.csv(climate_fossil,"climate_fossil.csv")

modern_pollen<- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv", row.names=1)
pred_core <- read_csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/pred_core.csv")
pred_t_core <- read_csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/pred_t_core.csv")
pred_core <- pred_core[,-1]
pred_t_core <- pred_t_core[,-1]

setwd("C:/Users/ml4418/Desktop/Master Project/Data/Plots/Climate changes whole")

plotdata<-pred_t_core
modern_gdd<-h1$pred_t_gdd;modern_alpha<-h1$pred_t_alpha;modern_Tmin<-h1$pred_t_Tmin

Plot_fossil(pred_t_core,h1$pred_t_gdd,h1$pred_t_alpha,h1$pred_t_Tmin,insol_summer_winter)

#Each site with standard error
pred_t_core$Site<-as.factor(pred_t_core$Site)
pred_t_core$Site<-fct_inorder(pred_t_core$Site, ordered = NA)
bci<-1.96
for(j in 1:nlevels(pred_t_core$Site)){
  sitename<-as.character(levels(pred_t_core$Site)[j])
  plotdata<-pred_t_core[which(pred_t_core$Site==sitename),]
  insol_site<-insol_summer_winter[which(insol_summer_winter$Site==sitename),]
  Plot_fossil_each_site_with_SE(plotdata,insol_site,sitename,bci)
}

# Relationships between reconstructed climates and climate drivers -------------------------------------
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Plots/Relationships")
plotdata<-pred_t_core;traindata<-pred_t
plot_relationship(plotdata)
plot_relationship1(plotdata,traindata)

  
# Gradient ----------------------------------------------------------------
pred_t_core <- read_csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/pred_t_core.csv")
pred_t_core <- pred_t_core[,-1]
plotdata<-pred_t_core[-which(pred_t_core$Site=="Zonar"|pred_t_core$Site=="Villarquemado"),1:6]
h11<-h1[-which(h1$name=="Zonar"|h1$name=="Villarquemado"),]

gradient_env_t<-rbind.data.frame(get_gradient(plotdata,500,500),get_gradient(plotdata,1500,500),get_gradient(plotdata,2500,500),get_gradient(plotdata,3500,500),get_gradient(plotdata,4500,500),get_gradient(plotdata,5500,500),get_gradient(plotdata,6500,500),get_gradient(plotdata,7500,500),get_gradient(plotdata,8500,500),get_gradient(plotdata,9500,500),get_gradient(plotdata,10500,500),get_gradient(plotdata,11500,500))

gradient_env_t$longitude<-h11$longitude;gradient_env_t$latitude<-h11$latitude;gradient_env_t$altitude<-h11$altitude
gradient_env_t$mean_gdd<-as.numeric(gradient_env_t$mean_gdd)
gradient_env_t$mean_alpha<-as.numeric(gradient_env_t$mean_alpha)
gradient_env_t$mean_Tmin<-as.numeric(gradient_env_t$mean_Tmin)
gradient_env_t$age<-as.numeric(gradient_env_t$age)
gradient_env_t$age<-factor(gradient_env_t$age,labels=c("0.5 ka","1.5 ka","2.5 ka","3.5 ka","4.5 ka","5.5 ka","6.5 ka","7.5 ka","8.5 ka","9.5 ka","10.5 ka","11.5 ka"))

gradient_env_t$altitudelabel<-NA
for(j in 1:nrow(gradient_env_t)){
  if(as.numeric(gradient_env_t[j,"altitude"])>1800){
    gradient_env_t[j,"altitudelabel"]<-"high"
  }else{
    gradient_env_t[j,"altitudelabel"]<-"low"
  }
}
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Plots/Gradient")
plotGradient<-gradient_env_t;method<-"Tolerance weigthed WA-PLS"

plot_spatial_gradient_colour(gradient_env_t,"Tolerance weighted WA-PLS")




