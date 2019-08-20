#Prepare the dataset

# Prepare modern pollen ---------------------------------------------------

setwd("C:/Users/ml4418/Desktop/Master Project/Data/Input data")
d <- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Feb2019_a_pct.csv", row.names=1)
taxaColMin<-which(colnames(d)=="Abies")
taxaColMax<-which(colnames(d)=="Zygophyllaceae")

#MI transformation to alpha
MI<-d$MI;w<-3;fai<-1/MI;F<-1+fai-(1+fai^w)^(1/w)
alpha<-1.26*MI*F

hist(MI)
hist(alpha)
plot(alpha~MI)

#Diversity
Hill_N1<-rep(NA,nrow(d))
for(j in 1:nrow(d)){
  p<-t(d[j,taxaColMin:taxaColMax])
  lnp<-log(p)
  plnp<-p*lnp
  Hill_N1[j]<-exp(-sum(plnp,na.rm=TRUE))
}


modern_pollen<-cbind.data.frame(d[,1:taxaColMin-1],alpha,Hill_N1,d[,taxaColMin:taxaColMax])
write.csv(modern_pollen,"Modern_Pollen_gdd_alpha_Tmin.csv")


# Prepare fossil pollen ---------------------------------------------------
Holocene <- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Holocene.csv", row.names=1)
taxaColMin<-which(colnames(Holocene)=="Abies")
taxaColMax<-which(colnames(Holocene)=="Zygophyllaceae")

#Diversity
Hill_N1<-rep(NA,nrow(Holocene))
for(j in 1:nrow(Holocene)){
  p<-t(Holocene[j,3:ncol(Holocene)])
  lnp<-log(p)
  plnp<-p*lnp
  Hill_N1[j]<-exp(-sum(plnp,na.rm=TRUE))
}

Holocene<-cbind.data.frame(Holocene[,c(1,2)],Hill_N1,Holocene[,taxaColMin:taxaColMax])
write.csv(Holocene,"Holocene.csv")


# Get the modern climate at the core sites --------------------------------------------------------------
h1 <- read_excel("C:/Users/ml4418/Desktop/Master Project/Data/Holocene data_clean_amalgam_relative_abundance.xls.xlsx", sheet = 1)
modern_pollen<- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv", row.names=1)
pred <- read_csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/pred.csv")
pred_t <- read_csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/pred_t.csv")
pred <- pred[,-1]
pred_t <- pred_t[,-1]

h1<-h1[-which(h1$name=="Marbore"),]
h1[which(h1[,"name"]=="Estanya"),"name"]<-"Estanya Catena"

h1$entity.number<-c(5557,NA,496,5518,5561,5559,NA)

h1<-cbind.data.frame(h1,modern_pollen[as.numeric(h1$entity.number),1:9])
h1<-cbind.data.frame(h1,pred[as.numeric(h1$entity.number),c("pred_gdd","pred_alpha","pred_Tmin")],pred_t[as.numeric(h1$entity.number),c("pred_gdd","pred_alpha","pred_Tmin")])
colnames(h1)[c(ncol(h1)-2,ncol(h1)-1,ncol(h1))]<-c("pred_t_gdd","pred_t_alpha","pred_t_Tmin")

setwd("C:/Users/ml4418/Desktop/Master Project/Data/Input data")
write.csv(h1,"h1.csv")