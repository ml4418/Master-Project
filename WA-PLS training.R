#WA-PLS training
rm(list=ls())
install.packages("reshape2")
install.packages("SDMTools")
install.packages("rioja")
install.packages("ggplot2")
install.packages("dplyr")

source('C:/Users/ml4418/Desktop/Master Project/Script/Tolerance weighted WA-PLS/Define functions.R', encoding = 'UTF-8')
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Input data")

modern_pollen<- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv", row.names=1)

taxaColMin<-which(colnames(modern_pollen)=="Abies")
taxaColMax<-which(colnames(modern_pollen)=="Zygophyllaceae")

########################################################################################################
#WA-PLS
component<-3

pred<-WAPLS_training(modern_pollen,taxaColMin,taxaColMax,component) #Use the self defined function

pred$Entity.name<-modern_pollen$Entity.name
pred$Hill_N1<-modern_pollen$Hill_N1
write.csv(pred,"pred.csv")

########################################################################################################
#Tolerance weighted WA-PLS
component<-2
P<-read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Optimum and Tolerance of each taxa.csv", row.names=1)

pred_t<-Tolerance_weighted_WAPLS_training(modern_pollen,taxaColMin,taxaColMax,P,component) #Use the self defined function

pred_t$Entity.name<-modern_pollen$Entity.name
pred_t$Hill_N1<-modern_pollen$Hill_N1
write.csv(pred_t,"pred_t.csv")

# Standard error by DKL ---------------------------------------------------
train<-modern_pollen[,c(taxaColMin:taxaColMax)]

SE_DKL_train<-cbind.data.frame(Tolerance_weighted_WAPLS_SE_DKL(train,P$t1),Tolerance_weighted_WAPLS_SE_DKL(train,P$t2),Tolerance_weighted_WAPLS_SE_DKL(train,P$t3))
colnames(SE_DKL_train)<-c("SE_gdd_DKL_train","SE_alpha_DKL_train","SE_Tmin_DKL_train")
pred_t<-cbind.data.frame(pred_t,SE_DKL_train)

write.csv(SE_DKL_train,"SE_DKL_train.csv")
write.csv(pred_t,"pred_t.csv")

########################################################################################################
#Plot the training results
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Plots/Training results")

training_results<-pred;training_type<-"WA-PLS"
plot_training_results(training_results,training_type)

training_results<-pred_t;training_type<-"Tolerance weighted WA-PLS"
plot_training_results(training_results,training_type)


########################################################################################################
#Get the fitness of the training results
fitness_original<-Training_fitness(pred) #Use the self defined function
fitness_with_t<-Training_fitness(pred_t) #Use the self defined function
fitness<-fitness_original
fitness<-fitness_with_t

paste("gdd: y=",fitness["a_gdd",1],"x+",fitness["b_gdd",1])
paste("alpha: y=",fitness["a_alpha",1],"x+",fitness["b_alpha",1])
paste("Tmin: y=",fitness["a_Tmin",1],"x+",fitness["b_Tmin",1])


