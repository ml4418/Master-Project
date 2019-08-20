#########################################################################################
##############################    Define functions  #####################################
#########################################################################################

# Optimum and tolerance of the taxa ---------------------------------------
#Define a function which can calculate the optimum (abundance-weighted mean) and tolerance (abundance-weighted standard deviation)
taxa_optimum_tolerance<-function(modern_pollen,taxaColMin,taxaColMax){
  library(reshape)
  modern_pollen<-melt(as.data.frame(modern_pollen),id=c(1:taxaColMin-1),measured=c(taxaColMin:taxaColMax))
  colnames(modern_pollen)[colnames(modern_pollen)=="variable"] <- "taxa"
  colnames(modern_pollen)[colnames(modern_pollen)=="value"] <- "abundance"
  modern_pollen$taxaCode<-as.numeric(modern_pollen$taxa)
  taxaName<-levels(modern_pollen[,"taxa"])
  nTaxa<-taxaColMax-taxaColMin+1
  
  #Abundance-weighted mean and standard deviation
  P<-data.frame()
  for(i in 1:nTaxa){
    print(i)
    modern_pollen_sub<-subset(modern_pollen,taxaCode==i)
    
    pk<-modern_pollen_sub$abundance
    x1k<-modern_pollen_sub$gdd
    x2k<-modern_pollen_sub$alpha
    x3k<-modern_pollen_sub$Tmin
    
    #climate space
    min_gdd<-min(modern_pollen_sub[,"gdd"][modern_pollen_sub[,"abundance"]>0]);max_gdd<-max(modern_pollen_sub[,"gdd"][modern_pollen_sub[,"abundance"]>0])
    min_alpha<-min(modern_pollen_sub[,"alpha"][modern_pollen_sub[,"abundance"]>0]);max_alpha<-max(modern_pollen_sub[,"alpha"][modern_pollen_sub[,"abundance"]>0])
    min_Tmin<-min(modern_pollen_sub[,"Tmin"][modern_pollen_sub[,"abundance"]>0]);max_Tmin<-max(modern_pollen_sub[,"Tmin"][modern_pollen_sub[,"abundance"]>0])
    
    #optimum
    u1<-sum(pk*x1k)/sum(pk)
    u2<-sum(pk*x2k)/sum(pk)
    u3<-sum(pk*x3k)/sum(pk)
    
    #sd
    sd1<-sd(modern_pollen_sub$gdd)
    sd2<-sd(modern_pollen_sub$alpha)
    sd3<-sd(modern_pollen_sub$Tmin)
    
    #tolerance
    t1<-sqrt(sum(pk*(x1k-u1)^2)/(2*sum(pk)))
    t2<-sqrt(sum(pk*(x2k-u2)^2)/(2*sum(pk)))
    t3<-sqrt(sum(pk*(x3k-u3)^2)/(2*sum(pk)))
    
    #abundance
    mean_abundance<-mean(modern_pollen_sub[,"abundance"])
    max_abundance<-max(modern_pollen_sub[,"abundance"])
    
    #the number of entities with non-zero abundance
    nonzero<-sum(modern_pollen_sub$abundance!=0)
    
    p<-cbind.data.frame(i,taxaName[i],mean_abundance,max_abundance,u1,sd1,t1,u2,sd2,t2,u3,sd3,t3,min_gdd,max_gdd,min_alpha,max_alpha,min_Tmin,max_Tmin,nonzero)
    P<-rbind.data.frame(P,p)
  }
  return(P)
}

# WA-PLS training ---------------------------------------------------------
WAPLS_training<-function(modern_pollen,taxaColMin,taxaColMax,component){
  library(rioja)
  train<- modern_pollen[,taxaColMin:taxaColMax]
  gdd<- modern_pollen$gdd
  alpha<- modern_pollen$alpha
  Tmin<- modern_pollen$Tmin
  
  mod1<-WAPLS(train,gdd,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  mod2<-WAPLS(train,alpha,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  mod3<-WAPLS(train,Tmin,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  
  pred<-cbind.data.frame(gdd,as.data.frame(fitted(mod1)[,component]),alpha,as.data.frame(fitted(mod2)[,component]),Tmin,as.data.frame(fitted(mod3)[,component]))
  names(pred)<-c("gdd","pred_gdd","alpha","pred_alpha","Tmin","pred_Tmin")
  
  pred$difference_gdd<-pred$pred_gdd-pred$gdd
  pred$difference_alpha<-pred$pred_alpha-pred$alpha
  pred$difference_Tmin<-pred$pred_Tmin-pred$Tmin
  return(pred)
}

# Tolerance weighted WA-PLS training --------------------------------------
Tolerance_weighted_WAPLS_training<-function(modern_pollen,taxaColMin,taxaColMax,P,component){
  library(rioja)
  train<- modern_pollen[,taxaColMin:taxaColMax]
  gdd<- modern_pollen$gdd
  alpha<- modern_pollen$alpha
  Tmin<- modern_pollen$Tmin
  
  mod1<-WAPLS(train,gdd,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  mod2<-WAPLS(train,alpha,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  mod3<-WAPLS(train,Tmin,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  
  sd_gdd<-t(P$sd1)
  sd_alpha<-t(P$sd2)
  sd_Tmin<-t(P$sd3)
  t_gdd<-t(P$t1)
  t_alpha<-t(P$t2)
  t_Tmin<-t(P$t3)
  
  #Give the tolerance-weight according to the equation derived from DKL
  train_with_t_gdd<-train
  train_with_t_alpha<-train
  train_with_t_Tmin<-train
  for(i in 1:ncol(train)){
    train_with_t_gdd[,i]<-train[,i]*sd_gdd[i]^2/(t_gdd[i]^2)
    train_with_t_alpha[,i]<-train[,i]*sd_alpha[i]^2/(t_alpha[i]^2)
    train_with_t_Tmin[,i]<-train[,i]*sd_Tmin[i]^2/(t_Tmin[i]^2)
  }
  
  
  #Results of tolerance weighted WA-PLS
  mod11 <- predict(mod1, newdata = train_with_t_gdd, sse = TRUE, nboot = 100, match.data = TRUE, verbose = TRUE)
  
  mod22 <- predict(mod2, newdata = train_with_t_alpha, sse = TRUE, nboot = 100, match.data = TRUE, verbose = TRUE)
  
  mod33 <- predict(mod3, newdata = train_with_t_Tmin, sse = TRUE, nboot = 100, match.data = TRUE, verbose = TRUE)
  
  pred_t<-cbind.data.frame(gdd,as.data.frame(mod11$fit)[,component],alpha,as.data.frame(mod22$fit)[,component],Tmin,as.data.frame(mod33$fit)[,component])
  names(pred_t)<-c("gdd","pred_gdd","alpha","pred_alpha","Tmin","pred_Tmin")
  
  pred_t$difference_gdd<-pred_t$pred_gdd-pred_t$gdd
  pred_t$difference_alpha<-pred_t$pred_alpha-pred_t$alpha
  pred_t$difference_Tmin<-pred_t$pred_Tmin-pred_t$Tmin
  return(pred_t)
}

# Fitness of training -----------------------------------------------------
Training_fitness<-function(training_results){
  lm1<-lm(pred_gdd~gdd,data=training_results)
  lm2<-lm(pred_alpha~alpha,data=training_results)
  lm3<-lm(pred_Tmin~Tmin,data=training_results)
  
  coef1<-t(coef(summary(lm1))[,"Estimate"])
  coef2<-t(coef(summary(lm2))[,"Estimate"])
  coef3<-t(coef(summary(lm3))[,"Estimate"])
  
  r.squared1<-1-sum((training_results$pred_gdd - training_results$gdd)^2)/sum((training_results$pred_gdd-mean(training_results$pred_gdd))^2)
  r.squared2<-1-sum((training_results$pred_alpha - training_results$alpha)^2)/sum((training_results$pred_alpha-mean(training_results$pred_alpha))^2)
  r.squared3<-1-sum((training_results$pred_Tmin - training_results$Tmin)^2)/sum((training_results$pred_Tmin-mean(training_results$pred_Tmin))^2)
  
  N<-nrow(training_results)
  RMSE1<-sqrt(sum((training_results$pred_gdd - training_results$gdd)^2/N))
  RMSE2<-sqrt(sum((training_results$pred_alpha - training_results$alpha)^2/N))
  RMSE3<-sqrt(sum((training_results$pred_Tmin - training_results$Tmin)^2/N))
  
  MAE1<-sum(abs(training_results$pred_gdd - training_results$gdd))/N
  MAE2<-sum(abs(training_results$pred_alpha - training_results$alpha))/N
  MAE3<-sum(abs(training_results$pred_Tmin - training_results$Tmin))/N
  
  summary_fitness<-t(round(cbind.data.frame(coef1,coef2,coef3,r.squared1,r.squared2,r.squared3,RMSE1,RMSE2,RMSE3,MAE1,MAE2,MAE3),digits=3))
  row.names(summary_fitness)<-c("b_gdd","a_gdd","b_alpha","a_alpha","b_Tmin","a_Tmin","r.squared_gdd","r.squared_alpha","r.squared_Tmin","RMSE_gdd","RMSE_alpha","RMSE_Tmin","MAE_gdd","MAE_alpha","MAE_Tmin")
  return(summary_fitness)
}

# WAPLS core --------------------------------------------------------------
WAPLS_core<-function(modern_pollen,taxaColMin,taxaColMax,component,core){
  library(rioja)
  train<- modern_pollen[,taxaColMin:taxaColMax]
  gdd<- modern_pollen$gdd
  alpha<- modern_pollen$alpha
  Tmin<- modern_pollen$Tmin
  
  mod1<-WAPLS(train,gdd,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  mod2<-WAPLS(train,alpha,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  mod3<-WAPLS(train,Tmin,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  
  gdd <- predict(mod1, newdata = core, sse = TRUE, nboot = 100, match.data = TRUE, verbose = TRUE)
  alpha <- predict(mod2, newdata = core, sse = TRUE, nboot = 100, match.data = TRUE, verbose = TRUE)
  Tmin <- predict(mod3, newdata = core, sse = TRUE, nboot = 100, match.data = TRUE, verbose = TRUE)
  
  pred_core <-cbind.data.frame(as.data.frame(gdd$fit)[,component],as.data.frame(alpha$fit)[,component],as.data.frame(Tmin$fit)[,component],as.data.frame(gdd$SEP.boot)[,component],as.data.frame(alpha$SEP.boot)[,component],as.data.frame(Tmin$SEP.boot)[,component])
  names(pred_core)<-c("gdd","alpha","Tmin","SE_gdd","SE_alpha","SE_Tmin")
  
  return(pred_core)
}

# Tolerance weighted WA-PLS core ------------------------------------------
Tolerance_weighted_WAPLS_core<-function(modern_pollen,taxaColMin,taxaColMax,P,component,core){
  library(rioja)
  train<- modern_pollen[,taxaColMin:taxaColMax]
  gdd<- modern_pollen$gdd
  alpha<- modern_pollen$alpha
  Tmin<- modern_pollen$Tmin
  
  mod1<-WAPLS(train,gdd,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  mod2<-WAPLS(train,alpha,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  mod3<-WAPLS(train,Tmin,npls=3,iswapls=TRUE, standx=FALSE, lean=FALSE,check.data=TRUE)
  
  sd_gdd<-t(P$sd1)
  sd_alpha<-t(P$sd2)
  sd_Tmin<-t(P$sd3)
  t_gdd<-t(P$t1)
  t_alpha<-t(P$t2)
  t_Tmin<-t(P$t3)
  
  
  core_with_t_gdd<-core
  core_with_t_alpha<-core
  core_with_t_Tmin<-core
  
  for(i in 1:ncol(core)){
    core_with_t_gdd[,i]<-core[,i]*sd_gdd[i]^2/(t_gdd[i]^2)
    core_with_t_alpha[,i]<-core[,i]*sd_alpha[i]^2/(t_alpha[i]^2)
    core_with_t_Tmin[,i]<-core[,i]*sd_Tmin[i]^2/(t_Tmin[i]^2)
  }
  
  gdd <- predict(mod1, newdata = core_with_t_gdd, sse = TRUE, nboot = 100, match.data = TRUE, verbose = TRUE)
  
  alpha <- predict(mod2, newdata = core_with_t_alpha, sse = TRUE, nboot = 100, match.data = TRUE, verbose = TRUE)
  
  Tmin <- predict(mod3, newdata = core_with_t_Tmin, sse = TRUE, nboot = 100, match.data = TRUE, verbose = TRUE)
  
  pred_t_core <-cbind.data.frame(as.data.frame(gdd$fit)[,component],as.data.frame(alpha$fit)[,component],as.data.frame(Tmin$fit)[,component],as.data.frame(gdd$SEP.boot)[,component],as.data.frame(alpha$SEP.boot)[,component],as.data.frame(Tmin$SEP.boot)[,component])
  names(pred_t_core)<-c("gdd","alpha","Tmin","SE_gdd","SE_alpha","SE_Tmin")
  
  return(pred_t_core)
}

# Standard error ----------------------------------------------------------
WAPLS_SE_replace1000<-function(train, env, n=1000, component, core = core){
  # Make NA filled list of names for each taxon 
  predr<-rep(NA, nrow(core))
  library(rioja)
  # Make many sets, run WAPLS 
  replicate(n, {                                         # Do this n times...
    k<-sample(1:nrow(train), size=nrow(train), replace=TRUE) # Make list of row numbers by sampling with replacement
    train<-train[k,]                                         # Reorganise train obs in k order
    train<-train[,colSums(train)>0]                            # Strip out zero-sum cols
    mod<-WAPLS(train, env[k])                                # Apply WAPLS, with env also in k order
    pred<-as.data.frame(predict(mod,core))[,component]      # Make reconstruction
    predr<-pred
    predr
  })
}  
Tolerance_weighted_WAPLS_SE_replace1000<-function(train, env, sd_env, t_env,n=1000, component, core = core){
  # Make NA filled list of names for each taxon 
  predr<-rep(NA, nrow(core))
  library(rioja)
  # Make many sets, run WAPLS 
  replicate(n, {                                         # Do this n times...
    k<-sample(1:nrow(train), size=nrow(train), replace=TRUE) # Make list of row numbers by sampling with replacement
    train<-train[k,]                                         # Reorganise train obs in k order
    train<-train[,colSums(train)>0]                            # Strip out zero-sum cols
    mod<-WAPLS(train, env[k])                                # Apply WAPLS, with env also in k order
    core_with_t_env<-core
    for(i in 1:ncol(core)){
      core_with_t_env[,i]<-core[,i]*sd_env[i]^2/(t_env[i]^2)
    }
    pred<-as.data.frame(predict(mod,core_with_t_env))[,component]      # Make reconstruction
    predr<-pred
    predr
  })
}  
Tolerance_weighted_WAPLS_SE_DKL<-function(core,t_env){
  SE_env<-rep(NA,nrow(core))
  for(j in 1:nrow(core)){
    p<-core[j,]
    SE_env[j]<-sqrt(1/sum(p/t_env^2))
  }
  return(SE_env)
}

# Get the gradient --------------------------------------------------------
get_gradient<-function(plotdata,age,range){
  library(forcats)
  plotdata$Site<-as.factor(plotdata$Site)
  plotdata$Site<-fct_inorder(plotdata$Site, ordered = NA)
  mean_env<-data.frame(matrix(NA,ncol=8,nrow=nlevels(plotdata$Site)))
  for(m in 1:nlevels(plotdata$Site)){
    plotdata_each<-plotdata[which(plotdata$Site==levels(plotdata$Site)[m]),]
    
    mean_gdd<-mean(as.matrix(plotdata_each[which(plotdata_each$Age.cal.BP>age-range & plotdata_each$Age.cal.BP<age+range),"gdd"]),na.rm = TRUE)
    mean_alpha<-mean(as.matrix(plotdata_each[which(plotdata_each$Age.cal.BP>age-range & plotdata_each$Age.cal.BP<age+range),"alpha"]),na.rm = TRUE)
    mean_Tmin<-mean(as.matrix(plotdata_each[which(plotdata_each$Age.cal.BP>age-range & plotdata_each$Age.cal.BP<age+range),"Tmin"]),na.rm = TRUE)
    SD_gdd<-sd(as.matrix(plotdata_each[which(plotdata_each$Age.cal.BP>age-range & plotdata_each$Age.cal.BP<age+range),"gdd"]),na.rm = TRUE)
    SD_alpha<-sd(as.matrix(plotdata_each[which(plotdata_each$Age.cal.BP>age-range & plotdata_each$Age.cal.BP<age+range),"alpha"]),na.rm = TRUE)
    SD_Tmin<-sd(as.matrix(plotdata_each[which(plotdata_each$Age.cal.BP>age-range & plotdata_each$Age.cal.BP<age+range),"Tmin"]),na.rm = TRUE)
    
    mean_env[m,]<-c(age,levels(plotdata$Site)[m],mean_gdd,mean_alpha,mean_Tmin,SD_gdd,SD_alpha,SD_Tmin)
  }
  colnames(mean_env)<-c("age","site","mean_gdd","mean_alpha","mean_Tmin","SD_gdd","SD_alpha","SD_Tmin")
  return(mean_env)
}

# Plot training results ---------------------------------------------------
plot_training_results<-function(training_results,training_type){
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(gridExtra)
  
  p1<-ggplot(training_results,aes(Tmin,pred_Tmin))+geom_point(size=0.8)+ labs(subtitle = "(a)")+
    geom_abline(slope=1, intercept=0)+labs(y= expression(paste("Reconstructed T"[min] ( degree~C))), x = expression(paste("T"[min] ( degree~C))))
  
  p2<-ggplot(training_results,aes(gdd,pred_gdd))+geom_point(size=0.8)+labs(subtitle = "(b)")+
    geom_abline(slope=1, intercept=0)+labs(y= bquote('Reconstructed'~ GDD[0]), x = bquote(GDD[0]))
  
  p3<-ggplot(training_results,aes(alpha,pred_alpha))+geom_point(size=0.8)+ labs(subtitle = "(c)")+
    geom_abline(slope=1, intercept=0)+labs(y= expression("Reconstructed "*alpha), x = expression(alpha))
  
  p4<-ggplot(training_results,aes(Tmin,difference_Tmin))+geom_point(size=0.8)+ labs(subtitle = "(d)")+
    geom_abline(slope=0, intercept=0)+labs(y= expression(paste("Residual of T"[min] ( degree~C))), x = expression(paste("T"[min] ( degree~C))))
  
  p5<-ggplot(training_results,aes(gdd,difference_gdd))+geom_point(size=0.8)+ labs(subtitle = "(e)")+
    geom_abline(slope=0, intercept=0)+labs(y= bquote('Residual of'~ GDD[0]), x = bquote(GDD[0]))
  
  p6<-ggplot(training_results,aes(alpha,difference_alpha))+geom_point(size=0.8)+ labs(subtitle = "(f)")+
    geom_abline(slope=0, intercept=0)+labs(y= expression("Residual of "*alpha), x = expression(alpha))
  
  
  p<-arrangeGrob(p1,p4,p2,p5,p3,p6 ,ncol = 2,top=paste(training_type,"training results"))
  ggsave(file=paste(training_type,"training results.png"),p,width=8,height=8)
}

# Plot the climate changes during the Holocene ----------------------------
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

Plot_fossil<-function(plotdata,modern_gdd,modern_alpha,modern_Tmin,insol_summer_winter){
  library(ggplot2)
  library(RColorBrewer)
  library(forcats)
  library(grid)
  library(gridExtra)
  
  plotdata$Site<-as.factor(plotdata$Site)
  plotdata$Site<-fct_inorder(plotdata$Site, ordered = NA)
  
  xbreak<-c(0,2000,4000,6000,8000,10000,12000)
  p1<-ggplot(plotdata, aes(x = Age.cal.BP, y = Tmin, colour = Site))+geom_point(size=0.8)+geom_line()+
    geom_vline(xintercept=8200)+geom_vline(xintercept=4200)+xlab("Age (cal yr BP)")+ylab(expression(paste("T"[min] (degree~C))))+
    labs(subtitle = "(a)")+scale_color_brewer(palette="Dark2")+
    scale_x_continuous(breaks = xbreak,limits = c(-1200,12000))+
    geom_segment(x = 0,xend=(-1000),y = as.numeric(modern_Tmin[1]),yend = as.numeric(modern_Tmin[1]), colour = brewer.pal(n = 7, name = "Dark2")[1])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_Tmin[2]),yend = as.numeric(modern_Tmin[2]), colour = brewer.pal(n = 7, name = "Dark2")[2])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_Tmin[3]),yend = as.numeric(modern_Tmin[3]),colour = brewer.pal(n = 7, name = "Dark2")[3])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_Tmin[5]),yend = as.numeric(modern_Tmin[5]), colour = brewer.pal(n = 7, name = "Dark2")[5])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_Tmin[6]),yend = as.numeric(modern_Tmin[6]), colour = brewer.pal(n = 7, name = "Dark2")[6])
  
  p2<-ggplot(plotdata, aes(x = Age.cal.BP, y = gdd, colour = Site))+geom_point(size=0.8)+geom_line()+
    geom_vline(xintercept=8200)+geom_vline(xintercept=4200)+xlab("Age (cal yr BP)")+ylab(bquote(GDD[0]))+
    labs(subtitle = "(b)")+ scale_color_brewer(palette="Dark2")+
    scale_x_continuous(breaks = xbreak,limits = c(-1200,12000))+
    geom_segment(x = 0,xend=(-1000),y = as.numeric(modern_gdd[1]),yend = as.numeric(modern_gdd[1]), colour = brewer.pal(n = 7, name = "Dark2")[1])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_gdd[2]),yend = as.numeric(modern_gdd[2]), colour = brewer.pal(n = 7, name = "Dark2")[2])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_gdd[3]),yend = as.numeric(modern_gdd[3]),colour = brewer.pal(n = 7, name = "Dark2")[3])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_gdd[5]),yend = as.numeric(modern_gdd[5]), colour = brewer.pal(n = 7, name = "Dark2")[5])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_gdd[6]),yend = as.numeric(modern_gdd[6]), colour = brewer.pal(n = 7, name = "Dark2")[6])
  
  p3<-ggplot(plotdata, aes(x = Age.cal.BP, y = alpha, colour = Site))+geom_point(size=0.8)+geom_line()+
    geom_vline(xintercept=8200)+geom_vline(xintercept=4200)+xlab("Age (cal yr BP)")+ylab(expression(alpha))+
    labs(subtitle = "(c)")+scale_color_brewer(palette="Dark2")+
    scale_x_continuous(breaks = xbreak,limits = c(-1200,12000))+ylim(0,1.05)+
    geom_segment(x = 0,xend=(-1000),y = as.numeric(modern_alpha[1]),yend = as.numeric(modern_alpha[1]), colour = brewer.pal(n = 7, name = "Dark2")[1])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_alpha[2]),yend = as.numeric(modern_alpha[2]), colour = brewer.pal(n = 7, name = "Dark2")[2])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_alpha[3]),yend = as.numeric(modern_alpha[3]),colour = brewer.pal(n = 7, name = "Dark2")[3])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_alpha[5]),yend = as.numeric(modern_alpha[5]), colour = brewer.pal(n = 7, name = "Dark2")[5])+
    geom_segment(x = 0,xend=(-1000), y = as.numeric(modern_alpha[6]),yend = as.numeric(modern_alpha[6]), colour = brewer.pal(n = 7, name = "Dark2")[6])
  
  
  insol_summer_winter$Site<-as.factor(insol_summer_winter$Site)
  insol_summer_winter$Site<-fct_inorder(insol_summer_winter$Site, ordered = NA)
  
  p4<-ggplot(insol_summer_winter)+geom_line(aes(Age,Insolation_winter,colour=Site))+
    labs(subtitle = "(d)")+ scale_color_brewer(palette="Dark2")+
    scale_x_continuous(breaks = xbreak,limits = c(0,12000))+
    ylab(bquote('Winter insolation (W'~m^-2~')'))
  
  p5<-ggplot(insol_summer_winter)+geom_line(aes(Age,Insolation_summer,colour=Site))+
    labs(subtitle = "(e)")+ scale_color_brewer(palette="Dark2")+
    scale_x_continuous(breaks = xbreak,limits = c(0,12000))+
    ylab(bquote('Summer insolation (W'~m^-2~')'))
  
  lengend<-get_legend(p1)
  p1<-p1+ theme(legend.position="none")
  p2<-p2+ theme(legend.position="none")
  p3<-p3+ theme(legend.position="none")
  p4<-p4+ theme(legend.position="none")
  p5<-p5+ theme(legend.position="none")
  
  p<-grid.arrange(p1,p2,p3,p4,p5,legend,layout_matrix=cbind(c(1,2,3),c(4,6,5)),
                  top="Climate reconstruction at the fossil sites")
  ggsave(file="Climate reconstruction at the fossil sites.png",p,width=9,height=10)
}

Plot_fossil_each_site_with_SE<-function(plotdata,insol_site,sitename,bci){
  library(ggplot2)
  library(RColorBrewer)
  library(forcats)
  library(grid)
  library(gridExtra)
  
  xbreak<-c(0,2000,4000,6000,8000,10000,12000)
  p1<-ggplot(plotdata, aes(x = Age.cal.BP, y = Tmin))+
    geom_errorbar(aes(ymin=Tmin-bci*SE_Tmin_DKL, ymax=Tmin+bci*SE_Tmin_DKL), colour="grey",width=0.2)+
    geom_point(aes(x = Age.cal.BP, y = Tmin),size=0.8)+geom_line()+
    geom_vline(xintercept=8200)+geom_vline(xintercept=4200)+xlab("Age (cal yr BP)")+ylab(expression(paste("T"[min] (degree~C))))+
    labs(subtitle = "(a)")+
    scale_x_continuous(breaks = xbreak,limits = c(-1200,12000))
  
  p2<-ggplot(plotdata, aes(x = Age.cal.BP, y = gdd))+
    geom_errorbar(aes(ymin=gdd-bci*SE_gdd_DKL, ymax=gdd+bci*SE_gdd_DKL), colour="grey", width=0.2)+
    geom_point(aes(x = Age.cal.BP, y = gdd),size=0.8)+geom_line()+
    geom_vline(xintercept=8200)+geom_vline(xintercept=4200)+xlab("Age (cal yr BP)")+ylab(bquote(GDD[0]))+
    labs(subtitle = "(b)")+ 
    scale_x_continuous(breaks = xbreak,limits = c(-1200,12000))
  
  p3<-ggplot(plotdata, aes(x = Age.cal.BP, y = alpha))+
    geom_errorbar(aes(ymin=alpha-bci*SE_alpha_DKL, ymax=alpha+bci*SE_alpha_DKL), colour="grey", width=0.2)+ 
    geom_point(aes(x = Age.cal.BP, y = alpha),size=0.8)+geom_line()+
    geom_vline(xintercept=8200)+geom_vline(xintercept=4200)+xlab("Age (cal yr BP)")+ylab(expression(alpha))+
    labs(subtitle = "(c)")+ 
    scale_x_continuous(breaks = xbreak,limits = c(-1200,12000))
  
  p4<-ggplot(insol_site)+geom_line(aes(Age,Insolation_summer))+
    geom_line(aes(Age,Insolation_winter*2.4),linetype="dashed")+
    labs(subtitle = "(d)")+
    scale_x_continuous(breaks = xbreak,limits = c(0,12000))+
    ylab(bquote('Summer insolation (W'~m^-2~')'))+
    scale_y_continuous(sec.axis = sec_axis(~.*(1/2.4), name = bquote('Winter insolation (W'~m^-2~')')))
  
  p<-arrangeGrob(p1,p2,p3,p4,ncol = 2,nrow=2,top=paste("Climate reconstruction at",sitename))
  ggsave(file=paste("Climate reconstruction at",sitename,".png"),p,width=10,height=6)
}

# Plot the relationships between reconstructed climates and climate drivers -------------------------------------
plot_relationship<-function(plotdata,traindata){
  library(ggplot2)
  library(RColorBrewer)
  library(forcats)
  library(grid)
  library(gridExtra)
  
  plotdata$Site<-as.factor(plotdata$Site)
  plotdata$Site<-fct_inorder(plotdata$Site, ordered = NA)
  
  p1<-ggplot(traindata, aes(x = alpha, y = gdd))+geom_point(size=0.8)+labs(subtitle = "(a)")+
    xlab(expression(alpha))+ylab(bquote(GDD[0]))
  
  p2<-ggplot(plotdata, aes(x = alpha, y = gdd,colour=Site))+geom_point(size=0.8)+labs(subtitle = "(b)")+
    xlab(expression(alpha))+ylab(bquote(GDD[0]))+ scale_color_brewer(palette="Dark2")
  
  p<-grid.arrange(p1,p2,ncol=2,widths=c(2,3),top="Relationships between warmth and dryness")
  ggsave(file="Relationships between gdd and alpha.png",p,width=8,height=3)
}

# Plot spatial gradient  ------------------------------------------
plot_spatial_gradient_colour<-function(plotGradient,method){
  
  plotGradient$altitudelabel <- factor(plotGradient$altitudelabel, levels = c("low","high"))
  p1<-ggplot(plotGradient)+labs(subtitle = expression(paste("(a) T"[min] ( degree~C))))+
    geom_point(aes(longitude,latitude,colour=mean_Tmin,shape=altitudelabel),size=3)+facet_grid(age~.)+
    scale_color_gradientn(colours = c("purple","blue","green","yellow","orange","red"),na.value = "grey90")+
    xlim(-5,5)+labs(colour = expression("Tmin" ( degree~C)),shape="altitude")+ theme(legend.position="bottom")+guides(shape = FALSE)+
    scale_y_continuous(breaks = c(42,43),limits = c(41,44))
  
  p2<-ggplot(plotGradient)+labs(subtitle = bquote('(b)'~ GDD[0]))+
    geom_point(aes(longitude,latitude,colour=mean_gdd,shape=altitudelabel),size=3)+facet_grid(age~.)+
    scale_color_gradientn(colours = c("purple","blue","green","yellow","orange","red"),na.value = "grey90")+
    xlim(-5,5)+labs(colour = bquote(GDD[0]))+ theme(legend.position="bottom")+guides(shape = FALSE)+
    scale_y_continuous(breaks = c(42,43),limits = c(41,44))+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank())
  
  p3<-ggplot(plotGradient)+labs(subtitle = expression("(c) "*alpha))+
    geom_point(aes(longitude,latitude,colour=mean_alpha,shape=altitudelabel),size=3)+facet_grid(age~.)+
    scale_color_gradientn(colours = c("yellow","green","blue","purple"),na.value = "grey90")+
    xlim(-5,5)+labs(colour = expression(alpha),shape="altitude")+ theme(legend.position="bottom")+guides(shape = FALSE)+
    scale_y_continuous(breaks = c(42,43),limits = c(41,44))+ 
    theme(axis.title.y = element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank()) 
  
  
  p<-arrangeGrob(p1,p2,p3,ncol = 3,top="Changes of gradient during the Holocene")
  ggsave(file="Changes of gradient during the Holocene.png",p,width=7,height=9)
}



########################################################################################################
############################    Prepare the dataset needed (Table 1) ###################################
########################################################################################################
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

########################################################################################################
#########################    Calculate optimum and tolerance (Table 3) #################################
########################################################################################################
# Calculate optimum and tolerance --------
library(reshape)
library(SDMTools)
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Input data")
modern_pollen<- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv", row.names=1)
taxaColMin<-which(colnames(modern_pollen)=="Abies")
taxaColMax<-which(colnames(modern_pollen)=="Zygophyllaceae")

P<-taxa_optimum_tolerance(modern_pollen,taxaColMin,taxaColMax)#Use the self defined function

write.csv(P,"Optimum and Tolerance of each taxa.csv")

# Non-zero samples (Table 3) --------------------------------------------------------
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




########################################################################################################
########################    PLots to show the principle (Figure 1)    ##################################
########################################################################################################
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Plots/Training results")
#Plots to show the principle
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

#u1>u2, t1>t2
a1_a<-log(0.6);u1_a<-10;t1_a<-8
a2_a<-log(0.6);u2_a<-(-6);t2_a<-4

x<-2;P1_a<-exp(a1_a-(x-u1_a)^2/(2*t1_a^2));P2_a<-exp(a2_a-(x-u2_a)^2/(2*t2_a^2))

pred_x<-(P1_a*u1_a+P2_a*u2_a)/(P1_a+P2_a);pred_xt<-(P1_a*u1_a/(t1_a^2)+P2_a*u2_a/(t2_a^2))/(P1_a/(t1_a^2)+P2_a/(t2_a^2))
SE<-sqrt(1/sum(P1_a/(t1_a^2)+P2_a/(t2_a^2)))

p1<-ggplot(data.frame(x=c(-30,30)), aes(x)) + 
  stat_function(fun=function(x) exp(a1_a-(x-u1_a)^2/(2*t1_a^2)))+
  stat_function(fun=function(x) exp(a2_a-(x-u2_a)^2/(2*t2_a^2)))+
  geom_vline(xintercept=pred_x,linetype="dashed")+geom_vline(xintercept=pred_xt,linetype="dashed")+
  annotate("text", x=u1_a,y=0,label = deparse(bquote(u[1])),parse = T)+
  annotate("text", x=u2_a,y=0,label = deparse(bquote(u[2])),parse = T)+
  annotate("text", x=pred_x,y=0.55,label = deparse(bquote(hat(x)[WA])),parse = T)+
  annotate("text", x=pred_xt,y=0.55,label = deparse(bquote(hat(x))),parse = T)+
  labs(subtitle = expression(paste("(a) ", u[1],">",u[2]," , ",t[1],">",t[2])))+theme_void()

#u1<u2, t1<t2
a2_b<-log(0.6);u2_b<-10;t2_b<-8
a1_b<-log(0.6);u1_b<-(-6);t1_b<-4

x<-2;P1_b<-exp(a1_b-(x-u1_b)^2/(2*t1_b^2));P2_b<-exp(a2_b-(x-u2_b)^2/(2*t2_b^2))

pred_x<-(P1_b*u1_b+P2_b*u2_b)/(P1_b+P2_b);pred_xt<-(P1_b*u1_b/(t1_b^2)+P2_b*u2_b/(t2_b^2))/(P1_b/(t1_b^2)+P2_b/(t2_b^2))
SE<-sqrt(1/sum(P1_b/(t1_b^2)+P2_b/(t2_b^2)))

p2<-ggplot(data.frame(x=c(-30,30)), aes(x)) + 
  stat_function(fun=function(x) exp(a1_b-(x-u1_b)^2/(2*t1_b^2)))+
  stat_function(fun=function(x) exp(a2_b-(x-u2_b)^2/(2*t2_b^2)))+
  geom_vline(xintercept=pred_x,linetype="dashed")+geom_vline(xintercept=pred_xt,linetype="dashed")+
  annotate("text", x=u1_b,y=0,label = deparse(bquote(u[1])),parse = T)+
  annotate("text", x=u2_b,y=0,label = deparse(bquote(u[2])),parse = T)+
  annotate("text", x=pred_x,y=0.55,label = deparse(bquote(hat(x)[WA])),parse = T)+
  annotate("text", x=pred_xt,y=0.55,label = deparse(bquote(hat(x))),parse = T)+
  labs(subtitle = expression(paste("(b) ", u[1],"<",u[2]," , ",t[1],"<",t[2])))+theme_void()

#u1>u2, t1<t2
a1_d<-log(0.6);u1_d<-10;t1_d<-4
a2_d<-log(0.6);u2_d<-(-6);t2_d<-8

x<-4;P1_d<-exp(a1_d-(x-u1_d)^2/(2*t1_d^2));P2_d<-exp(a2_d-(x-u2_d)^2/(2*t2_d^2))

pred_x<-(P1_d*u1_d+P2_d*u2_d)/(P1_d+P2_d);pred_xt<-(P1_d*u1_d/(t1_d^2)+P2_d*u2_d/(t2_d^2))/(P1_d/(t1_d^2)+P2_d/(t2_d^2))
SE<-sqrt(1/sum(P1_d/(t1_d^2)+P2_d/(t2_d^2)))

p3<-ggplot(data.frame(x=c(-30,30)), aes(x)) + 
  stat_function(fun=function(x) exp(a1_d-(x-u1_d)^2/(2*t1_d^2)))+
  stat_function(fun=function(x) exp(a2_d-(x-u2_d)^2/(2*t2_d^2)))+
  geom_vline(xintercept=pred_x,linetype="dashed")+geom_vline(xintercept=pred_xt,linetype="dashed")+
  annotate("text", x=u1_d,y=0,label = deparse(bquote(u[1])),parse = T)+
  annotate("text", x=u2_d,y=0,label = deparse(bquote(u[2])),parse = T)+
  annotate("text", x=pred_x,y=0.55,label = deparse(bquote(hat(x)[WA])),parse = T)+
  annotate("text", x=pred_xt,y=0.55,label = deparse(bquote(hat(x))),parse = T)+
  labs(subtitle = expression(paste("(c) ", u[1],">",u[2]," , ",t[1],"<",t[2])))+theme_void()


#u1<u2, t1>t2
a2_c<-log(0.6);u2_c<-10;t2_c<-4
a1_c<-log(0.6);u1_c<-(-6);t1_c<-8

x<-4;P1_c<-exp(a1_c-(x-u1_c)^2/(2*t1_c^2));P2_c<-exp(a2_c-(x-u2_c)^2/(2*t2_c^2))

pred_x<-(P1_c*u1_c+P2_c*u2_c)/(P1_c+P2_c);pred_xt<-(P1_c*u1_c/(t1_c^2)+P2_c*u2_c/(t2_c^2))/(P1_c/(t1_c^2)+P2_c/(t2_c^2))
SE<-sqrt(1/sum(P1_c/(t1_c^2)+P2_c/(t2_c^2)))

p4<-ggplot(data.frame(x=c(-30,30)), aes(x)) + 
  stat_function(fun=function(x) exp(a1_c-(x-u1_c)^2/(2*t1_c^2)))+
  stat_function(fun=function(x) exp(a2_c-(x-u2_c)^2/(2*t2_c^2)))+
  geom_vline(xintercept=pred_x,linetype="dashed")+geom_vline(xintercept=pred_xt,linetype="dashed")+
  annotate("text", x=u1_c,y=0,label = deparse(bquote(u[1])),parse = T)+
  annotate("text", x=u2_c,y=0,label = deparse(bquote(u[2])),parse = T)+
  annotate("text", x=pred_x,y=0.55,label = deparse(bquote(hat(x)[WA])),parse = T)+
  annotate("text", x=pred_xt,y=0.55,label = deparse(bquote(hat(x))),parse = T)+
  labs(subtitle = expression(paste("(d) ", u[1],"<",u[2]," , ",t[1],">",t[2])))+theme_void()

p_principle<-arrangeGrob(p1,p2,p3,p4,ncol = 2)

# Optimum and tolerance ---------------------------------------------------
P<-read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Optimum and Tolerance of each taxa.csv", row.names=1)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
p_Tmin<-ggplot(P,aes(u3,t3))+geom_point(size=0.8)+labs(subtitle = "(e)")+
  labs(y= expression(paste("Tolerance of T"[min] ( degree~C))), x =expression(paste("Optimum of T"[min] ( degree~C))))
p_gdd<-ggplot(P,aes(u1,t1))+geom_point(size=0.8)+labs(subtitle = "(f)")+
  labs(y= bquote('Tolerance of'~ GDD[0]), x = bquote('Optimum of'~ GDD[0]))
p_alpha<-ggplot(P,aes(u2,t2))+geom_point(size=0.8)+labs(subtitle = "(g)")+
  labs(y= expression("Tolerance of "*alpha), x = expression("Optimum of "*alpha))

p_optimum_tolerance<-arrangeGrob(p_Tmin,p_gdd,p_alpha,ncol = 1)

# Put them together -------------------------------------------------------
p<-arrangeGrob(p_principle,p_optimum_tolerance,ncol = 3,nrow=1,layout_matrix = cbind(1,1,2),padding=unit(2,"line"),top=textGrob("The principle of compression in WA-PLS", gp=gpar(fontsize=13)))

ggsave(file="The principle of compression in WA-PLS.png",p,width=9,height=6)



########################################################################################################
############    PLots to show the training data (Figure 2) and fossil data (Figure 3, Table 1)    ###############
########################################################################################################
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


########################################################################################################
#######################################    Training  ###################################################
########################################################################################################
# WA-PLS training ---------------------------------------------------------
rm(list=ls())
install.packages("reshape2")
install.packages("SDMTools")
install.packages("rioja")
install.packages("ggplot2")
install.packages("dplyr")

setwd("C:/Users/ml4418/Desktop/Master Project/Data/Input data")

modern_pollen<- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv", row.names=1)

taxaColMin<-which(colnames(modern_pollen)=="Abies")
taxaColMax<-which(colnames(modern_pollen)=="Zygophyllaceae")


# WA-PLS ------------------------------------------------------------------
component<-3

pred<-WAPLS_training(modern_pollen,taxaColMin,taxaColMax,component) #Use the self defined function

pred$Entity.name<-modern_pollen$Entity.name
pred$Hill_N1<-modern_pollen$Hill_N1
write.csv(pred,"pred.csv")


# Tolerance weighted WA-PLS -----------------------------------------------
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
#######################   Plots to show the training results (Figure 4, Table 2)  ###############################
########################################################################################################
# Plot the training results -----------------------------------------------
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Plots/Training results")
#Supporting information Figure 8
training_results<-pred;training_type<-"WA-PLS"
plot_training_results(training_results,training_type)

#Figure 4
training_results<-pred_t;training_type<-"Tolerance weighted WA-PLS"
plot_training_results(training_results,training_type)  

# Get the fitness of the training results (Table 2)---------------------------------
fitness_original<-Training_fitness(pred) #Use the self defined function
fitness_with_t<-Training_fitness(pred_t) #Use the self defined function
fitness<-fitness_original
fitness<-fitness_with_t

paste("gdd: y=",fitness["a_gdd",1],"x+",fitness["b_gdd",1])
paste("alpha: y=",fitness["a_alpha",1],"x+",fitness["b_alpha",1])
paste("Tmin: y=",fitness["a_Tmin",1],"x+",fitness["b_Tmin",1])


########################################################################################################
####################################    Reconstruction  ################################################
########################################################################################################
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



########################################################################################################
###########################    Calculate insolation at each site  #######################################
########################################################################################################
# Insolation script --------------------------------------------------------------
# purled from test_palinsol_01.Rmd
install.packages("palinsol")
library(palinsol)

# generate times
time_BP <- seq(-31000,28,by=1) #Where I made changes
head(time_BP); tail(time_BP)

# orbital parameters using Berger (1978) solution
orbital_params <- data.frame(time_BP, t(sapply(time_BP, function(tt) astro(tt, ber78, degree=TRUE))))
orbital_params$climprec <- orbital_params$ecc*sin(orbital_params$varpi*(2*pi/360.0))
orbital_params$perih360dyr <- as.numeric(t(sapply(time_BP, function (tt) 
  date_of_perihelion(astro(tt, ber78, degree=FALSE)))))
head(orbital_params); tail(orbital_params)

# plot the orbital parameters
plot(orbital_params$time_BP, orbital_params$eps, type="o", pch=16, main="Obliquity")
plot(orbital_params$time_BP, orbital_params$ecc, type="o", pch=16, main="Eccentricty")
plot(orbital_params$time_BP, orbital_params$climprec, type="o", pch=16, main="Climatic Precession")
plot(orbital_params$time_BP, orbital_params$varpi, type="o", pch=16, main="True Solar Longitude of Perihelion")
plot(orbital_params$time_BP, orbital_params$perih360dyr, type="o", pch=16, main="Day of Perihelion (in 360-Day Year)")
outfile="orbital_params_150kyr_1kyr.csv"
names(orbital_params) <- c("time_BP", "Obliq_degrees", "Eccentricity", "Lon_Perih_degrees", "epsp", 
                           "ClimaticPrecession", "Perihelion_360dyr")
write.table(orbital_params, outfile, sep=",", row.names=FALSE)

# Milankovitch plots -- present-day day of year vs. latitude
tt = 0.0
orbit <- astro(tt, ber78, degree = FALSE)
M <- Milankovitch(orbit)
str(M)

# plot using month as x-axis
plot(M, plot=contour, main=paste("Milankovitch Plot (t =", toString(tt), ")", sep=""))

# plot using true solar longitude as x-axis
plot(M, plot=contour, month=FALSE, main=paste("Milankovitch Plot (t =", toString(tt), ")", sep=""))

# Milankovitch plot of insolation long-term differences ("anomalies")
# present day insolation values
tt_present = 0.0
orbit_present <- astro(tt_present, ber78, degree = FALSE)
M_present <- Milankovitch(orbit_present) 

# paleo (6 ka) insolation values
tt_paleo = -6000.0
orbit_paleo <- astro(tt_paleo, ber78, degree = FALSE)
M_paleo <- Milankovitch(orbit_paleo)

# insolation differences
M <- M_paleo - M_present
plot(M, plot=contour, main=paste("Insolation Long-term Mean Difference (Paleo - Present) t =", toString(tt_paleo), sep=" "))

library(RColorBrewer)

# extract attributes
tsl <- attr(M, "long")
Col  = c(which(tsl >= (360-80)) , which(tsl < (360-80)))
Col = c(Col, Col[1])
MM = M[ Col,]
Month = c(tsl, tsl[1]+360) # day numbers
Latitude <- attr(M, "lat")

insol_diff <- as.array(matrix(t(MM[,]), nrow = length(Latitude), ncol=length(Month)))

cutpts <- c(-50,-40,-30,-20,-10,0,10,20,30,40,50)

image(Month, Latitude, t(insol_diff), axes=FALSE, breaks=cutpts, col=rev(brewer.pal(10,"PuOr")),
      main=paste("Insolation Long-term Mean Difference (Paleo - Present) t =", toString(tt_paleo), sep=" "))
contour(Month, Latitude, t(insol_diff), axes=FALSE, add=TRUE)
axis (1, at=seq(0,11)*30+15, labels=c('J','F','M','A','M','J','J','A','S','O','N','D'), tick=TRUE)
axis(2, at=seq(-90,90,30), labels=c('90S','60S', '30S','Eq.','30N','60N','90N'))
axis (3, at=seq(0,11)*30+15, labels=rep('',12))
axis(4, at=seq(-90,90,30), labels=rep('',7))

# generate present-day mid-month day numbers
mid_month_days_0ka <- seq(15.5, 345.5, by=30)
mid_month_days_0ka

# get true solar longitudes (TSL) for present-day mid-month day numbers
tt_present = 0.0
orbit_present <- astro(tt_present, ber78, degree = FALSE)
mid_month_tsl_0ka <- day2l(orbit_present, mid_month_days_0ka) #/(pi/180)
mid_month_tsl_0ka

# get orbital parameters for past 31 ka
time_BP <- seq(-31000,28,by=1)  #Where I made changes
orbital_params <- data.frame(time_BP, t(sapply(time_BP, function(tt) astro(tt, ber78, degree=FALSE))))
head(orbital_params)
tail(orbital_params)

# define output variables
length(time_BP)
insol_month <- matrix(0, nrow=length(time_BP), ncol=12)
insol_ann <- rep(0, length(time_BP))
insol_calins <- rep(0, length(time_BP))



# Insolation at core sites ------------------------------------------------
h1 <- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/h1.csv")
h1<-h1[,-1]
h1<-h1[order(h1$longitude),]

setwd("C:/Users/ml4418/Desktop/Master Project/Data/Insolation/Core sites insolation")
for(j in 1:nrow(h1)){
  Lat<-h1[j,"latitude"]
  # mid-month insolation values
  for (month in seq(1:12)) {
    tsl <- mid_month_tsl_0ka[month] 
    insol_month[,month] <- as.numeric(Insol(orbital_params, long=tsl, lat=Lat*pi/180, S0=1365))
  }
  
  # plot monthly insolation
  month <- 7
  plot(insol_month[,month] ~ time_BP, pch=16, type="o", main=paste("Insolation month ",toString(month), 
                                                                   ", Latitude = ", toString(Lat), sep=""))
  
  # annual average insolation
  for (n in seq(1:length(time_BP))) {
    orbit <- astro(time_BP[n], ber78, degree=FALSE)
    insol_ann[n] <- as.numeric(Insol_l1l2(orbit, l1=0,l2=2*pi,lat=Lat*pi/180, avg=TRUE,ell=TRUE))
  }
  
  # plot annual average insolation
  plot(insol_ann ~ time_BP, pch=16, type="o", main=paste("Annual Insolation, Latitude = ", toString(Lat), sep=""))
  
  outfile=paste("insol_",j,Lat,"_31kyr_1yr.csv")
  insol_out <- data.frame(time_BP, insol_month, insol_ann)
  names(insol_out) <- c("time_BP", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Ann")
  write.table(insol_out, outfile, sep=",", row.names=FALSE)
}
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Insolation/Core sites insolation")

insol_1<-read.csv("insol_ 1 43.27196 _31kyr_1yr.csv")
insol_2<-read.csv("insol_ 2 37.48273 _31kyr_1yr.csv")
insol_3<-read.csv("insol_ 3 40.49386 _31kyr_1yr.csv")
insol_4<-read.csv("insol_ 4 42.79899 _31kyr_1yr.csv")
insol_5<-read.csv("insol_ 5 42.54527 _31kyr_1yr.csv")
insol_6<-read.csv("insol_ 6 42.02826 _31kyr_1yr.csv")
insol_7<-read.csv("insol_ 7 42.12478 _31kyr_1yr.csv")

#Winter isnolation
insol_1_winter<-as.data.frame(rowMeans(insol_1[,c("Dec","Jan","Feb")]))
insol_2_winter<-as.data.frame(rowMeans(insol_2[,c("Dec","Jan","Feb")]))
insol_3_winter<-as.data.frame(rowMeans(insol_3[,c("Dec","Jan","Feb")]))
insol_4_winter<-as.data.frame(rowMeans(insol_4[,c("Dec","Jan","Feb")]))
insol_5_winter<-as.data.frame(rowMeans(insol_5[,c("Dec","Jan","Feb")]))
insol_6_winter<-as.data.frame(rowMeans(insol_6[,c("Dec","Jan","Feb")]))
insol_7_winter<-as.data.frame(rowMeans(insol_7[,c("Dec","Jan","Feb")]))

insol_winter<-cbind.data.frame(-(insol_1[,c("time_BP")]),insol_1_winter,insol_2_winter,insol_3_winter,insol_4_winter,insol_5_winter,insol_6_winter,insol_7_winter)
names(insol_winter)<-c("Age",as.character(h1$name))
library(reshape)
insol_winter<-melt(as.data.frame(insol_winter),id=1,measured=2:8)
names(insol_winter)<-c("Age","Site","Insolation_winter")

write.csv(insol_winter,"insol_winter.csv")

#Summer insolation
insol_1_summer<-as.data.frame(rowMeans(insol_1[,c("Jun","Jul","Aug")]))
insol_2_summer<-as.data.frame(rowMeans(insol_2[,c("Jun","Jul","Aug")]))
insol_3_summer<-as.data.frame(rowMeans(insol_3[,c("Jun","Jul","Aug")]))
insol_4_summer<-as.data.frame(rowMeans(insol_4[,c("Jun","Jul","Aug")]))
insol_5_summer<-as.data.frame(rowMeans(insol_5[,c("Jun","Jul","Aug")]))
insol_6_summer<-as.data.frame(rowMeans(insol_6[,c("Jun","Jul","Aug")]))
insol_7_summer<-as.data.frame(rowMeans(insol_7[,c("Jun","Jul","Aug")]))

insol_summer<-cbind.data.frame(-(insol_1[,c("time_BP")]),insol_1_summer,insol_2_summer,insol_3_summer,insol_4_summer,insol_5_summer,insol_6_summer,insol_7_summer)
names(insol_summer)<-c("Age",as.character(h1$name))
library(reshape)
insol_summer<-melt(as.data.frame(insol_summer),id=1,measured=2:8)
names(insol_summer)<-c("Age","Site","Insolation_summer")

write.csv(insol_summer,"insol_summer.csv")

#Put summer and winter insolation together
insol_summer_winter<-cbind.data.frame(insol_summer,insol_winter)
insol_summer_winter<-insol_summer_winter[,-c(4,5)]
write.csv(insol_summer_winter,"insol_summer_winter.csv")


# Put insolation data to the reconstruction results --------------------------------------------------------------
insol_site<-data.frame(matrix(nrow=nrow(pred_t_core),ncol=4))
for(j in 1:nrow(pred_t_core)){
  insol_site[j,1:4]<-insol_summer_winter[which(insol_summer_winter$Site==as.character(pred_t_core[j,"Site"])&insol_summer_winter$Age==as.numeric(pred_t_core[j,"Age.cal.BP"])),]
}
pred_t_core$insol_summer<-insol_site[,3]
pred_t_core$insol_winter<-insol_site[,4]
write.csv(pred_t_core,"pred_t_core.csv")



########################################################################################################
###########################    Show the reconstruction results  ########################################
########################################################################################################
# Plot the climates changes during the Holocene ---------------------------
#Get the modern climate of the fossil sites
library(readxl)
library(readr)
h1 <- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/h1.csv")
h1<-h1[,-1]
h1<-h1[order(h1$longitude),]

#Output for thesis (Table 4)
climate_fossil<-h1[,c("name","gdd","pred_t_gdd","alpha","pred_t_alpha","Tmin","pred_t_Tmin")]
climate_fossil[,c("gdd","pred_t_gdd")]<-round(climate_fossil[,c("gdd","pred_t_gdd")],digits=0)
climate_fossil[,c("alpha","pred_t_alpha")]<-round(climate_fossil[,c("alpha","pred_t_alpha")],digits=3)
climate_fossil[,c("Tmin","pred_t_Tmin")]<-round(climate_fossil[,c("Tmin","pred_t_Tmin")],digits=2)
write.csv(climate_fossil,"climate_fossil.csv")

#Plot the reconstruction results (Figure 5)
modern_pollen<- read.csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv", row.names=1)
pred_core <- read_csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/pred_core.csv")
pred_t_core <- read_csv("C:/Users/ml4418/Desktop/Master Project/Data/Input data/pred_t_core.csv")
pred_core <- pred_core[,-1]
pred_t_core <- pred_t_core[,-1]

setwd("C:/Users/ml4418/Desktop/Master Project/Data/Plots/Climate changes whole")
plotdata<-pred_t_core
modern_gdd<-h1$pred_t_gdd;modern_alpha<-h1$pred_t_alpha;modern_Tmin<-h1$pred_t_Tmin

Plot_fossil(pred_t_core,h1$pred_t_gdd,h1$pred_t_alpha,h1$pred_t_Tmin,insol_summer_winter)

#Each site with standard error (Supporting information Figure 9~15)
pred_t_core$Site<-as.factor(pred_t_core$Site)
pred_t_core$Site<-fct_inorder(pred_t_core$Site, ordered = NA)
bci<-1.96
for(j in 1:nlevels(pred_t_core$Site)){
  sitename<-as.character(levels(pred_t_core$Site)[j])
  plotdata<-pred_t_core[which(pred_t_core$Site==sitename),]
  insol_site<-insol_summer_winter[which(insol_summer_winter$Site==sitename),]
  Plot_fossil_each_site_with_SE(plotdata,insol_site,sitename,bci)
}

# Relationships between reconstructed climates and climate drivers (Figure 6)-------------------------------------
setwd("C:/Users/ml4418/Desktop/Master Project/Data/Plots/Relationships")
plotdata<-pred_t_core;traindata<-pred_t
plot_relationship(plotdata,traindata)


# Gradient (Figure 7)----------------------------------------------------------------
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


