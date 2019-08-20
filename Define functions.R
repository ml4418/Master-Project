
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
plot_relationship<-function(plotdata){
  library(ggplot2)
  library(RColorBrewer)
  library(forcats)
  library(grid)
  library(gridExtra)
  
  plotdata$Site<-as.factor(plotdata$Site)
  plotdata$Site<-fct_inorder(plotdata$Site, ordered = NA)
  
  p1<-ggplot(plotdata, aes(x = alpha, y = gdd,colour=Site))+geom_point(size=0.8)+labs(subtitle = "(a)")+
    xlab(expression(alpha))+ylab(bquote(GDD[0]))+ scale_color_brewer(palette="Dark2")
  lengend<-get_legend(p1)
  p1<-p1+ theme(legend.position="none")
  
  p2<-ggplot(plotdata, aes(x = alpha, y = Tmin,colour=Site))+geom_point(size=0.8)+labs(subtitle = "(b)")+
    xlab(expression(alpha))+ylab(paste("T"[min] ( degree~C))))+ scale_color_brewer(palette="Dark2")+
    theme(legend.position="none")
  
  p3<-ggplot(plotdata, aes(x = Tmin, y = gdd,colour=Site))+geom_point(size=0.8)+labs(subtitle = "(c)")+
    xlab(paste("T"[min] ( degree~C))))+ylab(bquote(GDD[0]))+ scale_color_brewer(palette="Dark2")+
    theme(legend.position="none")
  
  p<-grid.arrange(p1,p2,p3,legend, ncol=2,
                  top="Relationships between reconstructed climates")
  ggsave(file="Relationships between reconstructed climates.png",p,width=7,height=5)
}
plot_relationship1<-function(plotdata,traindata){
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



