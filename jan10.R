setwd("N:/Thesis")
data1=read.csv("all_trans_short_r.csv",header=TRUE) #load data with NA's
data1cut=read.csv("trans_short_cut_r.csv",header=TRUE)

#install packages
install.packages("nlme")
install.packages("lme4")
install.packages("dplyr")
install.packages("nortest")
install.packages("lattice")
install.packages("agricolae")
install.packages("chron")

library(nlme)
library(lme4)
library(dplyr)
library(nortest)
library(lattice)
library(agricolae)
library(chron)
source("AEDgraphingfunctions.R")

#site specific non seasonal variables
boxplot(data2$elev,      xlab="All Sites",ylab="Elevation (m)",cex.lab=1.5,pch=16,group=data2$basin.stream)
boxplot(data2$from.north,xlab="All Sites",ylab="Aspect (Deg. South Facing)",cex.lab=1.5,pch=16,group=data2$basin.stream)
boxplot(data2$slope,     xlab="All Sites",ylab="Slope (%)",cex.lab=1.5,pch=16,group=data2$basin.stream)
boxplot(data2$bf,        xlab="All Sites",ylab="Bank Full (m)",cex.lab=1.5,pch=16,group=data2$basin.stream)
boxplot(data2$pebble,    xlab="All Sites",ylab="Pebble Median (mm)",cex.lab=1.5,pch=16,group=data2$basin.stream)
# Variables that change over time
interaction.plot(data1$season.yr,data1$stream,data1$width,xlab="Sampling Period", ylab="Wetted Width (m)",cex.lab=1.5,col="black",lwd=2.5,legend=F)
plot(HSD.test(aov(data1$width~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Width")


interaction.plot(data1$season.yr,data1$stream,data1$depth,xlab="Sampling Period", ylab="Depth (m)",cex.lab=1.5,col="black",lwd=2.5,legend=F)
plot(HSD.test(aov(data1$depth~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Depth")

interaction.plot(data1$season.yr,data1$stream,data1$discharge,xlab="Sampling Period", ylab="Discharge (L/s)",cex.lab=1.5,col="black",lwd=2.5,legend=F)
plot(HSD.test(aov(data1$discharge~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Discharge")
summary(aov(data1$discharge~data1$stream+data1$season.yr))

plot(data1$depth~data1$width,xlab = "Wetted Width (m)", ylab = "Depth (m)",cex.lab=1.5,pch=16)
abline(lm(data1$depth~data1$width))
summary(lm(data1$depth~data1$width))
legend("topleft", legend=c("R2= 0.69", "P= 9.1e-9"),bty="n")

plot(data1$discharge~data1$width,xlab = "Wetted Width (m)", ylab = "Discharge (L/s)",cex.lab=1.5,pch=16)
abline(lm(data1$discharge~data1$width))
summary(lm(data1$discharge~data1$width))
legend("topleft", legend=c("R2= 0.47", "P= 2.0e-5"),bty="n")

plot(data1$cut.mass~data1$cut.pop.dens)


interaction.plot(data1$season.yr,data1$stream,data1$carbon,xlab="Sampling Period", ylab="Carbon (DOC mg C/L)",ylim=c(0,14),cex.lab=1.5,col="black",lwd=2.5,legend=F)
legend("topleft", legend=c("P= 1.2e-4"),bty="n")
plot(HSD.test(aov(data1$carbon~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Carbon")
anova(aov(data1$carbon~data1$stream+data1$season.yr))


interaction.plot(data1$season.yr,data1$stream,data1$phosphate,xlab="Sampling Period", ylab="Phosphate (SRP mg P/L)",ylim=c(0,.06),cex.lab=1.5,col="black",lwd=2.5,legend=F)
legend("topleft", legend=c("P= 1.7e-7"),bty="n")
plot(HSD.test(aov(data1$phosphate~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Phosphate")
anova(aov(data1$phosphate~data1$stream+data1$season.yr))


interaction.plot(data1$season.yr,data1$stream,data1$din.out,xlab="Sampling Period", ylab="Nitrogen (DIN mg N/L)*",ylim=c(0,.017),cex.lab=1.5,col="black",lwd=2.5,legend=F) # added.016478 to ammonia and .0012613 to nitrate, dropped outliers
legend("bottomright", legend=c("P= 4.6e-7"),bty="n")
plot(HSD.test(aov(data1$din.out~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="DIN")
anova(aov(data1$din.out~data1$stream+data1$season.yr))

interaction.plot(data1$season.yr,data1$stream,data1$canopy,xlab="Sampling Period", ylab="Canopy (% open)",ylim=c(0,80),cex.lab=1.5,col="black",lwd=2.5,legend=F)
plot(HSD.test(aov(data1$canopy~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Canopy")
summary(aov(data1$canopy~data1$stream+data1$season.yr))


interaction.plot(data1$season.yr,data1$stream,data1$par.integrative,xlab="Sampling Period", ylab="PAR (integrative mol/m2/d)",ylim=c(0,3.5),cex.lab=1.5,col="black",lwd=2.5,legend=F) #integrative seemed more informative than mean or max PAR
legend("topright", legend=c("P= 0.0010"),bty="n")
plot(HSD.test(aov(data1$par.integrative~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="PAR")
anova(aov(data1$par.integrative~data1$stream+data1$season.yr))


interaction.plot(data1$season.yr,data1$stream,data1$temp.min,xlab="Sampling Period", ylab="Minimum Temp (deg C)",ylim=c(0,11),cex.lab=1.5,col="black",lwd=2.5,legend=F) # stream daily minimum temp seemed a better predictor than mean or max temp
legend("topright", legend=c("P= 2.6e-11"),bty="n")
plot(HSD.test(aov(data1$temp.min~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Min Temp")
anova(aov(data1$temp.min~data1$stream+data1$season.yr))

interaction.plot(data1cut$season.yr,data1cut$stream,data1cut$cut.mass.m,xlab="Sampling Period", ylab="Trout Biomass (g/m)",ylim=c(0,14),cex.lab=1.5,col="black",lwd=2.5,legend=F) # almost entirely CUTT with EBT in Jack summer '18 mixed in
plot(HSD.test(aov(data1cut$cut.mass.m~data1cut$stream+data1cut$season.yr), 'data1cut$season.yr'),ylab="Trout Biomass")
summary(aov(data1cut$cut.mass.m~data1cut$stream+data1cut$season.yr))

#metabolism modeling
# insert "example" graphs, may need to load things from "example" tab
plot(StreamData$dtime, StreamData$light,cex.lab=1.5,xaxt="n",ylim=c(0,50000),xlab="Time (Midnight to Midnight)",ylab="PAR (mol/m2/s)")
plot(StreamData2$dtime, StreamData2$sat,cex.lab=1.5,xaxt="n",ylim=c(91.75,94),xlab="Time (Midnight to Midnight)",ylab="O2 (% saturation)")
plot(StreamData$dtime, StreamData$temp,cex.lab=1.5,xaxt="n",ylim=c(7,11),xlab="Time (Midnight to Midnight)",ylab="Temp (Deg. C)")
plot(StreamData$dtime, StreamData$oxy,cex.lab=1.5,xaxt="n",ylim=c(9,9.8),xlab="Time (Midnight to Midnight)",ylab="O2 (mg/L)")
rivermetabK(o2file=StreamData[ StreamData $dtime>=as.numeric(chron(dates="07/01/18", times="00:00:00")) & StreamData $dtime<=as.numeric(chron(dates="07/01/18", times="23:59:59")), ], z=0.050, bp=710, ts=0.003472)
# metabolism.lit vs metabolism.my
plot(k600.my~k600.lit,xlab = "K600 Literature (1/d)", ylab = "K600 My Data",cex.lab=1.5,pch=16,data=data1)
abline(lm(k600.my~k600.lit,data=data1))
legend("topright", legend=c("R2= 0.058", "P= 0.11"),bty="n")
summary(lm(k600.my~k600.lit,data=data1))

plot(gpp.my~gpp.lit,xlab = "GPP Literature (g O2/m2/d)", ylab = "GPP My Data",cex.lab=1.5,pch=16,data=data1)
abline(lm(gpp.my~gpp.lit,data=data1))
legend("topleft", legend=c("R2= 0.54", "P= 9.1e-05"),bty="n")
summary(lm(gpp.my~gpp.lit,data=data1))

plot(abs(er.my)~abs(er.lit),xlab = "ER Literature (g O2/m2/d)", ylab = "ER My Data",cex.lab=1.5,pch=16,data=data1)
abline(lm(abs(er.my)~abs(er.lit),data=data1))
legend("topleft", legend=c("R2= 0.46", "P= 4e-4"),bty="n")
summary(lm(abs(er.my)~abs(er.lit),data=data1))

plot(pr.my~pr.lit,xlab = "PR ratio Literature", ylab = "PR Ratio My Data",cex.lab=1.5,pch=16,data=data1)
abline(lm(pr.my~pr.lit,data=data1))
legend("topleft", legend=c("R2= 0.80", "P= 3.2e-8"),bty="n")
summary(lm(pr.my~pr.lit,data=data1))

plot(r2.my~r2.lit,xlab = "R2 of models Literature", ylab = "R2 of models My Data",cex.lab=1.5,pch=16,data=data1)
abline(lm(r2.my~r2.lit,data=data1))
legend("topleft", legend=c("R2= 0.96", "P= 2.1e-15"),bty="n")
summary(lm(r2.my~r2.lit,data=data1))

# metab.lit
boxplot(data1$gpp.lit~data1$season.yr,xlab="Sampling Period", ylab="GPP (g O2/m2/d)",cex.lab=1.5,pch=16)
legend("topright", legend=c("P= 0.020"),bty="n")
plot(HSD.test(aov(data1$gpp.lit~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="GPP.lit")
anova(aov(data1$gpp.lit~data1$stream+data1$season.yr))

boxplot(abs(data1$er.lit)~data1$season.yr,xlab="Sampling Period", ylab="ER (g O2/m2/d)*",cex.lab=1.5,pch=16)
legend("topright", legend=c("P= 0.052"),bty="n")
plot(HSD.test(aov(abs(data1$er.lit)~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="ER.lit")
anova(aov(abs(data1$er.lit)~data1$stream+data1$season.yr))

boxplot(data1$pr.lit~data1$season.yr,xlab="Sampling Period", ylab="PR Ratio",cex.lab=1.5,pch=16)
plot(HSD.test(aov(data1$pr.lit~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="PR.lit")

# pairplot most   
Z=cbind(data1$from.north,data1$slope,data1$bf,data1$pebble,data1$velocity.mean,data1$discharge,data1$width,data1$depth,data1$carbon,data1$din.out,data1$phosphate,data1$canopy,data1$par.integrative,data1$temp.max,data1$temp.min)
colnames(Z)<-c("Aspect",      "Slope",   "BF",  "Pebble",          "Velocity",           "Q",       "Width",    "Depth",        "C",         "N",            "P",     "Canopy",        "Light",            "T.max",     "T.min")
pairs(Z[,1:15], lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)

# pairplot for models  
Z=cbind(data1$basin,data1$stream,data1$season.yr,data1$slope,data1$width,data1$depth,data1$canopy,data1$temp.min,data1$t.gpp.lit,data1$t.er.lit,data1$t.cut.mass.m)
colnames(Z)<-c("Basin",  "Stream",   "seas.yr",      "Slope",     "Width",    "Depth",   "Canopy",      "Min Temp",       "t.GPP"  ,"t.ER",         "t.Fish")
pairs(Z[,1:11], lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)





########################### GPP
data1$t.gpp.lit=log(1+data1$gpp.lit)^(1/2)
vf1=varIdent(form = ~1|stream) # alternate variance structure
M1.1<-lme(t.gpp.lit~season.yr+depth, na.action=na.omit, random = ~1|stream, data=data1, weights=vf1)
x.gpp<-data1$t.gpp.lit[!is.na(data1$t.gpp.lit)]#removes na values from column
E1.1.gpp<-residuals(M1.1,type="normalized")
summary(M1.1)

qqnorm(E1.1.gpp,main="",ylab="GPP Quantiles",cex.lab=1.5,pch=16)
qqline(E1.1.gpp)
legend("topleft", legend=c("AD normality test", "P= 0.14"),bty="n")
ad.test(E1.1.gpp)

plot(M1.1,cex.lab=1.5,pch=16,col=25,xlab="GPP model fitted values",ylab="Standardized Residuals") 

plot(x.gpp, E1.1.gpp,xlab="GPP transformed",ylab="Normalized Residuals",cex.lab=1.5,pch=16)
abline(0,0)
legend("topleft", legend=c("R2= 0.20", "P= 0.013"),bty="n")
summary(lm(x.gpp~E1.1.gpp))

plot(t.gpp.lit~season.yr,xlab="Sampling Period",ylab="GPP Transformed",cex.lab=1.5,pch=16,data=data1)
legend("topright", legend=c("Model P< 0.0001"),bty="n")
anova(M1.1)

plot(t.gpp.lit~depth,xlab="Depth (m)",ylab="GPP Transformed",ylim=c(0,.75),xlim=c(.035,.13),cex.lab=1.5,pch=16,data=data1)
legend("topleft", legend=c("R2= 0.13", "Model P< 0.0001"),bty="n")
abline(lm(t.gpp.lit~depth,data=data1))
anova(M1.1)
#GPP non model graphs
plot(t.gpp.lit~carbon,xlab="C (DOC mg/L)",ylab="GPP Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.gpp.lit~carbon,data=data1))
legend("topright", legend=c("R2= -0.021", "P= 0.50"),bty="n")
summary(lm(t.gpp.lit~carbon,data=data1))

plot(t.gpp.lit~din.out,xlab="N (DIN mg/L)*",ylab="GPP Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.gpp.lit~din.out,data=data1))
legend("topright", legend=c("R2= 0.13", "P= 0.050"),bty="n")
summary(lm(t.gpp.lit~din.out,data=data1))

plot(t.gpp.lit~phosphate,xlab="P (SRP mg/L)",ylab="GPP Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.gpp.lit~phosphate,data=data1))
legend("topright", legend=c("R2= 0.033", "P= 0.19"),bty="n")
summary(lm(t.gpp.lit~phosphate,data=data1))

plot(t.gpp.lit~par.integrative,xlab="PAR (mol photons/m2/d)",ylab="GPP Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.gpp.lit~par.integrative,data=data1))
legend("topright", legend=c("R2= 0.012", "P= 0.26"),bty="n")
summary(lm(t.gpp.lit~par.integrative,data=data1))

############################### ER
data1$t.er.lit=log(1+abs(data1$er.lit))
Mer.2=lme(t.er.lit~depth+slope, na.action=na.omit, random = ~1|stream, data=data1, method="ML") #adds random effect
x.er<-data1$t.er.lit[!is.na(data1$t.er.lit)]#removes na values from column
E.2.er<-residuals(Mer.2,type="normalized")

summary(Mer.2)
which(E.2.er>2)
E.2.er


qqnorm(E.2.er,main="",ylab="ER Quantiles",cex.lab=1.5,pch=16)
qqline(E.2.er)
legend("topleft", legend=c("AD normality test", "P= 0.16"),bty="n")
ad.test(E.2.er)

plot(Mer.2,cex.lab=1.5,pch=16,col=25,xlab="ER model fitted values",ylab="Standardized Residuals") 

plot(x.er, E.2.er,xlab="ER transformed",ylab="Normalized Residuals",cex.lab=1.5,pch=16)
legend("topleft", legend=c("R2= 8.0e-4", "P= 0.32"),bty="n")
abline(0,0)
summary(lm(x.er~E.2.er))

plot(t.er.lit~depth,xlab="Depth (m)",ylab="ER Transformed",ylim=c(1.5,3.5),xlim=c(.04,.13),data=data1,cex.lab=1.5,pch=16)
abline(lm(t.er.lit~depth,data=data1))
legend("bottomright", legend=c("R2= 0.36", "Model P< 0.0001"),bty="n")
summary(lm(t.er.lit~depth,data=data1))
anova(Mer.2)

plot(t.er.lit~slope,xlab="Slope (%)",ylab="ER Transformed",ylim=c(1.5,3.5),xlim=c(1.5,10.5),data=data1,cex.lab=1.5,pch=16)
abline(lm(t.er.lit~slope,data=data1))
legend("topleft", legend=c("R2= 0.57", "Model P= <0.0001"),bty="n")
summary(lm(t.er.lit~slope,data=data1))
anova(Mer.2)

# ER non model plots
plot(t.er.lit~carbon,xlab="C (DOC mg/L)",ylab="ER Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.er.lit~carbon,data=data1))
legend("topright", legend=c("R2= -0.042", "P= 0.99"),bty="n")
summary(lm(t.er.lit~carbon,data=data1))

plot(t.er.lit~din.out,xlab="N (DIN mg/L)*",ylab="ER Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.er.lit~din.out,data=data1))
legend("topleft", legend=c("R2= 0.030", "P= 0.21"),bty="n")
summary(lm(t.er.lit~din.out,data=data1))

plot(t.er.lit~phosphate,xlab="P (SRP mg/L)",ylab="ER Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.er.lit~phosphate,data=data1))
legend("topright", legend=c("R2= -0.025", "P= 0.53"),bty="n")
summary(lm(t.er.lit~phosphate,data=data1))

plot(t.er.lit~par.integrative,xlab="PAR (mol/m2/d)",ylab="ER Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.er.lit~par.integrative,data=data1))
legend("topright", legend=c("R2= -0.030", "P= 0.61"),bty="n")
summary(lm(t.er.lit~par.integrative,data=data1))

plot(t.er.lit~t.gpp.lit,data=data1)
summary(lm(t.er.lit~t.gpp.lit,data=data1))

plot(abs(er.lit)~gpp.lit,data=data1)
summary(lm(abs(er.lit)~gpp.lit,data=data1))



############################### CUT
data1$t.cut.mass.m2=log(1+data1$cut.mass.m2)
x.cut<-data1$t.cut.mass.m2[!is.na(data1$cut.mass.m2)]#removes na values from column
vf3=varExp(form = ~fitted(.))
Mcut1.3<-gls(t.cut.mass.m2~temp.min*canopy+basin, na.action=na.omit, data=data1, weights=vf3)
summary(Mcut1.3)

qqnorm(E1.3.cut,main="",ylab="Trout Quantiles",cex.lab=1.5,pch=16)
legend("topleft", legend=c("AD normality test", "P= 0.22"),bty="n")
qqline(E1.3.cut)
ad.test(E1.3.cut)

plot(Mcut1.3,cex.lab=1.5,pch=16,col=25,xlab="Trout model fitted values",ylab="Standardized Residuals") # outlier? Hurley Summer '18

plot(x.cut, E1.3.cut,xlab="Trout Transformed",ylab="Normalized Residuals",cex.lab=1.5,pch=16)
abline(0,0)
legend("topleft", legend=c("R2= 0.065", "P= 0.14"),bty="n")
summary(lm(x.cut~E1.3.cut))


plot(t.cut.mass.m~basin,xlab="Basin",ylab="Trout Transformed",cex.lab=1.5,pch=16,data=data1)
legend("top", legend=c("Model P= 3e-4"),bty="n")
anova(Mcut1.3)

plot(t.cut.mass.m~temp.min,xlab="Temp. min (deg C)",ylab="Trout Transformed",xlim=c(6.5,11),ylim=c(0,3),cex.lab=1.5,pch=16,data=data1)
legend("topright", legend=c("R2= 0.19", "Model P= 0.0001"),bty="n")
abline(lm(t.cut.mass.m~temp.min,data=data1))
summary(lm(t.cut.mass.m~temp.min,data=data1))
anova(Mcut1.3)

plot(t.cut.mass.m~canopy,xlab="Canopy (% open)",ylab="Trout Transformed",xlim=c(0,60),ylim=c(0,3),cex.lab=1.5,pch=16,data=data1)
abline(lm(t.cut.mass.m~canopy,data=data1))
legend("bottomright", legend=c("R2= -0.033", "Model P= 0.026"),bty="n")
summary(lm(t.cut.mass.m~canopy,data=data1))
anova(Mcut1.3)

data1cut$t.cut.mass.m=log(1+data1cut$cut.mass.m)
interaction.plot(data1cut$temp.min.cat,data1cut$canopy.cat,data1cut$t.cut.mass.m,ylab="Trout Transformed",xlab="Temperature Category",cex.lab=1.5,col="black",lwd=2.5,legend=F ) 
legend("bottomleft", legend=c("Model P= 0.0018"),bty="n")
anova(Mcut1.3)


plot(t.cut.mass.m~width,xlab="Width (m)",ylab="Trout Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.cut.mass.m~width,data=data1))
legend("bottomright", legend=c("R2= 0.27", "Model P= 0.0007"),bty="n")
summary(lm(t.cut.mass.m~width,data=data1))
anova(Mcut1.3)

# CUT non model graphs
plot(t.cut.mass.m~carbon,xlab="C (DOC mg/L)",ylab="Trout Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.cut.mass.m~carbon,data=data1))
legend("topright", legend=c("R2= -0.027", "P= 0.49"),bty="n")
summary(lm(t.cut.mass.m~carbon,data=data1))

plot(t.cut.mass.m~din.out,xlab="N (DIN mg/L)*",ylab="Trout Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.cut.mass.m~din.out,data=data1))
legend("topleft", legend=c("R2= -0.051", "P= 0.74"),bty="n")
summary(lm(t.cut.mass.m~din.out,data=data1))

plot(t.cut.mass.m~phosphate,xlab="P (SRP mg/L)",ylab="Trout Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.cut.mass.m~phosphate,data=data1))
legend("bottomleft", legend=c("R2= 0.15", "P= 0.054"),bty="n")
summary(lm(t.cut.mass.m~phosphate,data=data1))

plot(t.cut.mass.m~par.integrative,xlab="Light (PAR mol/m2/d)",ylab="Trout Transformed",cex.lab=1.5,pch=16,data=data1)
abline(lm(t.cut.mass.m~par.integrative,data=data1))
legend("topright", legend=c("R2= -0.053", "P= 0.85"),bty="n")
summary(lm(t.cut.mass.m~par.integrative,data=data1))

### Other considerations
plot(data1$t.er.lit~data1$t.gpp.lit,xlab ="GPP Transformed" , ylab ="ER Transformed",cex.lab=1.5,pch=16)
abline(lm(data1$t.er.lit~data1$t.gpp.lit))
legend("topleft", legend=c("R2= 0.36", "P= 6.8e-4"),bty="n")
summary(lm(data1$t.er.lit~data1$t.gpp.lit))

plot(data1$er.lit~data1$gpp.lit,xlab ="GPP (g O2/m2/d)" , ylab ="ER (g O2/m2/d)",cex.lab=1.5,pch=16)
abline(lm(data1$er.lit~data1$gpp.lit))
legend("bottomleft", legend=c("R2= 0.41", "P= 2.6e-4"),bty="n")
summary(lm(data1$er.lit~data1$gpp.lit))

plot(data1$t.cut.mass.m~data1$t.gpp.lit,xlab ="GPP Transformed" , ylab ="Trout Transformed",cex.lab=1.5,pch=16)
abline(lm(data1$t.cut.mass.m~data1$t.gpp.lit))
legend("topright", legend=c("R2= -0.054", "P= 0.89"),bty="n")
summary(lm(data1$t.cut.mass.m~data1$t.gpp.lit))

plot(data1$t.cut.mass.m~data1$t.er.lit,xlab ="ER Transformed" , ylab ="Trout Transformed",cex.lab=1.5,pch=16)
abline(lm(data1$t.cut.mass.m~data1$t.er.lit))
summary(lm(data1$t.cut.mass.m~data1$t.er.lit))
legend("bottomright", legend=c("R2= -0.053", "P= 0.84"),bty="n")

plot(data1$cut.mass.m2~data1$gpp.lit,xlab ="GPP" , ylab ="Trout",cex.lab=1.5,pch=16)
abline(lm(data1$cut.mass.m2~data1$gpp.lit))
legend("topright", legend=c("R2= 0.012", "P= 0.28"),bty="n")
summary(lm(data1$cut.mass.m2~data1$gpp.lit))


plot(data1$carbon~data1$phosphate,xlab ="P (mg P/L)" , ylab ="C (DOC mg C/L)",cex.lab=1.5,pch=16)
abline(lm(data1$carbon~data1$phosphate))
summary(lm(data1$carbon~data1$phosphate))
legend("topright", legend=c("R2= 0.19", "P= 0.010"),bty="n")

x.carb=data1$carbon[!is.na(data1$carbon)]#removes na values from column
Mcarb=lm(data1$phosphate~data1$carbon)
Ecarb=residuals(Mcarb)

plot(data1$carbon,Ecarb,xlab="C (DOC mg C/L)",ylab="Residuals",cex.lab=1.5,pch=16)
abline(0,0)

