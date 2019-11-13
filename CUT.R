data1=read.csv("all_trans_short_r.csv",header=TRUE) #load data in short format with NA's
#install packages
install.packages("nlme")
install.packages("lme4")
install.packages("dplyr")
install.packages("nortest")


library(nlme)
library(lme4)
library(dplyr)
library(nortest)
source("AEDgraphingfunctions.R") #we used this alot in advanced biostats


# CUT mass/m transformation

dotchart(data1$cut.mass.m)
hist(data1$cut.mass.m)
ad.test(data1$cut.mass.m)
data1$t.cut.mass.m=log(1+data1$cut.mass.m)
ad.test(data1$t.cut.mass.m)
dotchart(data1$t.cut.mass.m)
hist(data1$t.cut.mass.m)

summary(lm(cut.mass.m~basin,na.action=na.omit,data=data1))          #p .02   r2 .3  ###
summary(lm(cut.mass.m~stream,na.action=na.omit,data=data1))         #p .2    r2 .2
summary(lm(cut.mass.m~yr,na.action=na.omit,data=data1))             #p .4    r2 -.004
summary(lm(cut.mass.m~elev,na.action=na.omit,data=data1))           #p .04   r2 .16 ###
summary(lm(cut.mass.m~from.north,na.action=na.omit,data=data1))     #p .8    r2 -.05
summary(lm(cut.mass.m~slope,na.action=na.omit,data=data1))          #p .2    r2 .04
summary(lm(cut.mass.m~bf,na.action=na.omit,data=data1))             #p .2    r2 .05
summary(lm(cut.mass.m~pebble,na.action=na.omit,data=data1))         #p .7    r2 -.05
summary(lm(cut.mass.m~width,na.action=na.omit,data=data1))          #p .02   r2 .2  ###
summary(lm(cut.mass.m~depth,na.action=na.omit,data=data1))          #p .18   r2 .05
summary(lm(cut.mass.m~velocity.mean,na.action=na.omit,data=data1))  #p .6    r2 -.03
summary(lm(cut.mass.m~discharge,na.action=na.omit,data=data1))      #p .3    r2 .02
summary(lm(cut.mass.m~carbon,na.action=na.omit,data=data1))         #p .8    r2 -.05
summary(lm(cut.mass.m~ammonia,na.action=na.omit,data=data1))        #p .3    r2 -.004
summary(lm(cut.mass.m~nitrate,na.action=na.omit,data=data1))        #p .2    r2 .04
summary(lm(cut.mass.m~din,na.action=na.omit,data=data1))            #p .2    r2 .04
summary(lm(cut.mass.m~phosphate,na.action=na.omit,data=data1))      #p .4    r2 -.01
summary(lm(cut.mass.m~cn.ratio,na.action=na.omit,data=data1))       #p .4    r2 -.01
summary(lm(cut.mass.m~canopy,na.action=na.omit,data=data1))         #p .3    r2 .001
summary(lm(cut.mass.m~par.integrative,na.action=na.omit,data=data1))#p .9    r2 -.05
summary(lm(cut.mass.m~par.mean,na.action=na.omit,data=data1))       #p .9    r2 -.05
summary(lm(cut.mass.m~par.max,na.action=na.omit,data=data1))        #p .16   r2 .06
summary(lm(cut.mass.m~temp.max,na.action=na.omit,data=data1))       #p .2    r2 .02
summary(lm(cut.mass.m~temp.min,na.action=na.omit,data=data1))       #p .04   r2 .16 ###
summary(lm(cut.mass.m~temp.mean,na.action=na.omit,data=data1))      #p .07   r2 .13 ###
summary(lm(cut.mass.m~gpp.lit,na.action=na.omit,data=data1))        #p .5    r2 -.03
summary(lm(cut.mass.m~pr.lit,na.action=na.omit,data=data1))         #p .5    r2 -.03
summary(lm(cut.mass.m~cut.capt,na.action=na.omit,data=data1))       #p .5    r2 -.03
summary(lm(cut.mass.m~cut.pop,na.action=na.omit,data=data1))        #p .004  r2 .3  ###
summary(lm(cut.mass.m~cut.mass,na.action=na.omit,data=data1))       #p .17   r2 .05
summary(lm(cut.mass.m~sculp.mass.m,na.action=na.omit,data=data1))   #p .08   r2 .11

plot(data1$cut.mass.m,data1$par.integrative)
ad.test(residuals(lmGPP))

Z=cbind(data1$basin,data1$width,data1$elev,data1$temp.min,data1$stream,data1$cut.mass.m)
colnames(Z)<-c("basin", "width",    "elev",  "temp.min",      "stream",  "cut biomass")
pairs(Z[,1:6], lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)

lmCUT1=lm(cut.mass.m~basin+width+elev+temp.min+stream,data=data1)
lmCUT1=lm(cut.mass.m~temp.min+elev+width,data=data1)

summary(lmCUT1)
lm(formula=lmCUT1)
anova(lmCUT1)
drop1(lmCUT1,test="F")
step(lmCUT1)

summary(lm(t.cut.mass.m~basin,na.action=na.omit,data=data1))          #p .004   r2 .4
summary(lm(t.cut.mass.m~stream,na.action=na.omit,data=data1))         #p .1    r2 .4
summary(lm(t.cut.mass.m~yr,na.action=na.omit,data=data1))             #p .5    r2 -.02
summary(lm(t.cut.mass.m~elev,na.action=na.omit,data=data1))           #p .3   r2 .2
summary(lm(t.cut.mass.m~from.north,na.action=na.omit,data=data1))     #p .8    r2 -.05
summary(lm(t.cut.mass.m~slope,na.action=na.omit,data=data1))          #p .18    r2 .05
summary(lm(t.cut.mass.m~bf,na.action=na.omit,data=data1))             #p .16    r2 .06
summary(lm(t.cut.mass.m~pebble,na.action=na.omit,data=data1))         #p .69    r2 -.05
summary(lm(t.cut.mass.m~width,na.action=na.omit,data=data1))          #p .013   r2 .3
summary(lm(t.cut.mass.m~depth,na.action=na.omit,data=data1))          #p .14   r2 .07
summary(lm(t.cut.mass.m~velocity.mean,na.action=na.omit,data=data1))  #p .48    r2 -.03
summary(lm(t.cut.mass.m~discharge,na.action=na.omit,data=data1))      #p .18    r2 .05
summary(lm(t.cut.mass.m~carbon,na.action=na.omit,data=data1))         #p .61    r2 -.04
summary(lm(t.cut.mass.m~ammonia,na.action=na.omit,data=data1))        #p .42    r2 -.02
summary(lm(t.cut.mass.m~nitrate,na.action=na.omit,data=data1))        #p .22    r2 .03
summary(lm(t.cut.mass.m~din,na.action=na.omit,data=data1))            #p .21    r2 .04
summary(lm(t.cut.mass.m~phosphate,na.action=na.omit,data=data1))      #p .18    r2 .05
summary(lm(t.cut.mass.m~cn.ratio,na.action=na.omit,data=data1))       #p .37    r2 -.008
summary(lm(t.cut.mass.m~canopy,na.action=na.omit,data=data1))         #p .43    r2 -.02
summary(lm(t.cut.mass.m~par.integrative,na.action=na.omit,data=data1))#p .86    r2 -.05
summary(lm(t.cut.mass.m~par.mean,na.action=na.omit,data=data1))       #p .86  r2 -.05
summary(lm(t.cut.mass.m~par.max,na.action=na.omit,data=data1))        #p .23   r2 .02
summary(lm(t.cut.mass.m~temp.max,na.action=na.omit,data=data1))       #p .22  r2 .03
summary(lm(t.cut.mass.m~temp.min,na.action=na.omit,data=data1))       #p .033   r2 .2
summary(lm(t.cut.mass.m~temp.mean,na.action=na.omit,data=data1))      #p .057   r2 .14
summary(lm(t.cut.mass.m~gpp.lit,na.action=na.omit,data=data1))        #p .62    r2 -.04
summary(lm(t.cut.mass.m~pr.lit,na.action=na.omit,data=data1))         #p .61  r2 -.04
summary(lm(t.cut.mass.m~cut.capt,na.action=na.omit,data=data1))       #p .72    r2 -.05
summary(lm(t.cut.mass.m~cut.pop,na.action=na.omit,data=data1))        #p .001  r2 .4
summary(lm(t.cut.mass.m~cut.mass,na.action=na.omit,data=data1))       #p .16   r2 .06
summary(lm(t.cut.mass.m~sculp.mass.m,na.action=na.omit,data=data1))   #p .13   r2 .07


Z=cbind(data1$basin,data1$width,data1$temp.min,data1$stream,data1$t.cut.mass.m)
colnames(Z)<-c("basin", "width",    "temp.min",  "stream",     "t.cut biomass")
pairs(Z[,1:5], lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)


lmCUT2=lm(t.cut.mass.m~basin+width+temp.min+stream,data=data1)

summary(lmCUT2)
qqnorm(E.1.cut)
qqline(E.1.cut)
ad.test(E.1.cut)
plot(Mcut.1) 
lmCUT2resid=resid(lmCUT2)
plot(x.cut,lmCUT2resid)

lm(formula=lmCUT2)
anova(lmCUT2)
drop1(lmCUT2,test="F")
step(lmCUT2)

#predictor transformations
dotchart(data1$temp.min)
data1$t.temp.min=log(1+data1$temp.min)
dotchart(data1$t.temp.min)

dotchart(data1$elev)
data1$t.elev=log(1+data1$elev)
dotchart(data1$t.elev)

dotchart(data1$width)
data1$t.width=log(1+data1$width)
dotchart(data1$t.width)

###############################################################################
#using glm to analyze variables
Mcut.1=gls(t.cut.mass.m~temp.min*canopy+width+basin,na.action=na.omit, 
           data=data1, method="ML") #base model
Mcut.2=lme(t.cut.mass.m~temp.min*canopy+width+basin, na.action=na.omit, 
           random = ~1|stream, data=data1, method="ML") #adds random effect
Mcut.3=lme(t.cut.mass.m~temp.min*canopy+width+basin, na.action=na.omit, 
           random = ~1|basin.stream, data=data1, method="ML") #adds nesting

anova(Mcut.1, Mcut.2, Mcut.3)

## # # Analyze residuals
x.cut<-data1$t.cut.mass.m[!is.na(data1$cut.mass.m)]#removes na values from column
# # Mcut.1 residuals
E.1.cut<-residuals(Mcut.1,type="normalized")
qqnorm(E.1.cut)
qqline(E.1.cut)
ad.test(E.1.cut)
plot(Mcut.1) 
plot(x.cut, E.1.cut)
summary(lm(x.cut~E.1.cut))

# # # # M.2 residuals
E.2.cut<-residuals(Mcut.2,type="normalized")
qqnorm(E.2.cut)
qqline(E.2.cut)
ad.test(E.2.cut)
plot(Mcut.2) 
plot(x.cut, E.2.cut)
summary(lm(x.cut~E.2.cut))
plot(data1$t.width[!is.na(data1$cut.mass.m)],E.2.cut)
plot(data1$t.temp.min[!is.na(data1$cut.mass.m)],E.2.cut)
plot(data1$basin[!is.na(data1$cut.mass.m)],E.2.cut)
abline(1,0)
abline(-1,0)
abline(lm(E.2.cut~x.cut))

# # # M.3 residuals
E.3.cut<-residuals(Mcut.3,type="normalized")
qqnorm(E.3.cut)
qqline(E.3.cut)
ad.test(E.3.cut)
plot(Mcut.3) 
plot(x.cut, E.3.cut)
summary(lm(x.cut~E.3.cut))

#residuals somewhat linear? try alternate variance structures?
#####################################################
#Analyze models with alternate variance structures
#####################################################
#base model from above
Mcut.1=gls(t.cut.mass.m~t.width+t.temp.min+basin, na.action=na.omit, 
           data=data1)#base model #run without "ML", method="REML" instead, default for comparison

#####alternate variance structures
vf1=varIdent(form = ~1|stream)
vf2=varPower(form = ~fitted(.))
vf3=varExp(form = ~fitted(.))
vf4=varConstPower(form = ~fitted(.))
####vf5=varPower(form = ~fitted (.)|site)
####vf6=varExp(form = ~fitted(.)|site)
####vf7=varConstPower(form = ~fitted(.)|site)
####vf8=varIdent(form = ~1|f.sample.event)
######alternate models
Mcut1.1<-gls(t.cut.mass.m~temp.min*canopy+width+basin, na.action=na.omit, 
             data=data1, weights=vf1)

Mcut1.2<-gls(t.cut.mass.m~temp.min*canopy+width+basin, na.action=na.omit, 
             data=data1, weights=vf2)

Mcut1.3<-gls(t.cut.mass.m~temp.min*canopy+width+basin, na.action=na.omit, 
             data=data1, weights=vf3)

Mcut1.4<-gls(t.cut.mass.m~temp.min*canopy+width+basin, na.action=na.omit, 
             data=data1, weights=vf4)

Mcut1.5<-gls(t.cut.mass.m~temp.min*canopy+width+basin, na.action=na.omit, 
             data=data1,  weights=varComb(vf1,vf2))

anova(Mcut.1,Mcut1.1,Mcut1.2,Mcut1.3,Mcut1.4)
anova(Mcut.1,Mcut1.1)
#Mcut1.1
#Analyze alternate variance structure model residuals
E1.1.cut<-residuals(Mcut1.1,type="normalized")
qqnorm(E1.1.cut)
qqline(E1.1.cut)
ad.test(E1.1.cut)
plot(Mcut1.1) 
plot(x.cut, E1.1.cut)
summary(lm(x.cut~E1.1.cut))
plot(E1.1.cut,data1cut$width)
plot(E1.1.cut,data1cut$temp.min)

E1.2.cut<-residuals(Mcut1.2,type="normalized")
qqnorm(E1.2.cut)
qqline(E1.2.cut)
ad.test(E1.2.cut)
plot(Mcut1.2) 
plot(x.cut, E1.2.cut)
summary(lm(x.cut~E1.2.cut))

E1.3.cut<-residuals(Mcut1.3,type="normalized")
qqnorm(E1.3.cut)
qqline(E1.3.cut)
ad.test(E1.3.cut)
plot(Mcut1.3) 
summary(Mcut1.3) 
plot(x.cut, E1.3.cut)
summary(lm(x.cut~E1.3.cut))

E1.4.cut<-residuals(Mcut1.4,type="normalized")
qqnorm(E1.4.cut)
qqline(E1.4.cut)
ad.test(E1.4.cut)
plot(Mcut1.4) 
plot(x.cut, E1.4.cut)
summary(lm(x.cut~E1.4.cut))

E1.5.cut<-residuals(Mcut1.5,type="normalized")
qqnorm(E1.5.cut)
qqline(E1.5.cut)
ad.test(E1.5.cut)
plot(Mcut1.5) 
plot(x.cut, E1.5.cut)
summary(lm(x.cut~E1.5.cut))

acf(E1.1.cut, na.action=na.pass,
    main="Auto-correlation plot for residuals")#check for autocorrelation in residuals


# CUT model validation
op=par(mfrow = c(2, 2)) 
plot(lmCUT2)
op=par(mfrow = c(2, 2)) 
CUT2resid=rstandard(lmCUT2) 
hist(CUT2resid) 
qqnorm(CUT2resid) 
plot(CUT2resid ~ data1cut$width, xlab = "width", ylab = "Residuals",cex.lab=1.5,pch=16) 
abline(0, 0)
plot(CUT2resid ~ data1cut$temp.min, xlab = "min temp", ylab = "Residuals",cex.lab=1.5,pch=16) 
abline(0, 0)


############################################## model interpretation
# CUT model interpretation
plot(data1cut$width,data1cut$t.cut.mass.m,cex.lab=1.5,pch=16)
abline(lm(data1cut$t.cut.mass.m~data1cut$width))
plot(data1cut$temp.min,data1cut$t.cut.mass.m,cex.lab=1.5,pch=16)
abline(lm(data1cut$t.cut.mass.m~data1cut$temp.min))
