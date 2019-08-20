#Calculate optimum (abundance-weighted mean) and tolerance (abundance-weighted standard deviation)
library(reshape)
library(SDMTools)
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Input data")
modern_pollen<- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv", row.names=1)
taxaColMin<-which(colnames(modern_pollen)=="Abies")
taxaColMax<-which(colnames(modern_pollen)=="Zygophyllaceae")

source('C:/Users/ml4418/Desktop/Master Project/Script/Tolerance weighted WA-PLS/Define functions.R', encoding = 'UTF-8')
P<-taxa_optimum_tolerance(modern_pollen,taxaColMin,taxaColMax)#Use the self defined function

write.csv(P,"Optimum and Tolerance of each taxa.csv")

# Non-zero samples --------------------------------------------------------
narrow_t_Tmin<-round(mean(P[which(P$t3<1.5),"nonzero"]),digits=0)
narrow_t_gdd<-round(mean(P[which(P$t1<250),"nonzero"]),digits=0)
narrow_t_alpha<-round(mean(P[which(P$t2<0.05),"nonzero"]),digits=0)

narrow_t_Tmin
narrow_t_gdd
narrow_t_alpha

nrow(modern_pollen)-narrow_t_Tmin
nrow(modern_pollen)-narrow_t_gdd
nrow(modern_pollen)-narrow_t_alpha

round(mean(P[,"nonzero"]),digits=0)
nrow(modern_pollen)-round(mean(P[,"nonzero"]),digits=0)
