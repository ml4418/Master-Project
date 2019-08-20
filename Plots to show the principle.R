# Principle code ----------------------------------------------------------
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
  