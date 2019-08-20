
# Insolation script I got --------------------------------------------------------------
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

ggplot(insol_winter,aes(Age,Insolation_winter,colour=Site))+geom_line()+
  scale_x_continuous(breaks = xbreak,limits = c(0,12000))+ scale_color_brewer(palette="Dark2")

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
ggplot(insol_summer,aes(Age,Insolation_summer,colour=Site))+geom_line()+
  scale_x_continuous(breaks = xbreak,limits = c(0,12000))+ scale_color_brewer(palette="Dark2")

#Put summer and winter insolation together
insol_summer_winter<-cbind.data.frame(insol_summer,insol_winter)
insol_summer_winter<-insol_summer_winter[,-c(4,5)]
write.csv(insol_summer_winter,"insol_summer_winter.csv")

ggplot(insol_summer_winter)+geom_line(aes(Age,Insolation_summer,colour=Site))+
  geom_line(aes(Age,Insolation_winter*2.4,colour=Site),linetype="dashed")+
  scale_x_continuous(breaks = xbreak,limits = c(0,12000))+
  scale_y_continuous(sec.axis = sec_axis(~.*(1/2.4), name = "Insolation_winter"))+ 
  scale_color_brewer(palette="Dark2")

