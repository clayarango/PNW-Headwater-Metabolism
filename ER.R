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


# ER transformation
dotchart(data1$er.lit,xlab="ER",cex.lab=1.5,pch=16)
hist(data1$er.lit,xlab="ER",cex.lab=1.5,pch=16)
ad.test(data1$er.lit)
data1$t.er.lit=log(1+abs(data1$er.lit))
ad.test(data1$t.er.lit)
hist(data1$t.er.lit,xlab="transfromed ER",cex.lab=1.5)
dotchart(data1$t.er.lit,xlab="transformed ER",cex.lab=1.5,pch=16)

# ER model selsection
summary(lm(t.er.lit~basin,na.action=na.omit,data=data1))          #p .01   r2 .2
summary(lm(t.er.lit~stream,na.action=na.omit,data=data1))         #p .9e-6 r2 .8
summary(lm(t.er.lit~season.yr,na.action=na.omit,data=data1))      #p .6    r2 -.04 
summary(lm(t.er.lit~yr,na.action=na.omit,data=data1))             #p .6    r2 -.03
summary(lm(t.er.lit~elev,na.action=na.omit,data=data1))           #p .4    r2 -.01
summary(lm(t.er.lit~from.north,na.action=na.omit,data=data1))     #p .9    r2 -.04
summary(lm(t.er.lit~slope,na.action=na.omit,data=data1))          #p .5e-6 r2 .6 
summary(lm(t.er.lit~bf,na.action=na.omit,data=data1))             #p .04   r2 .1
summary(lm(t.er.lit~pebble,na.action=na.omit,data=data1))         #p .001  r2 .3  
summary(lm(t.er.lit~width,na.action=na.omit,data=data1))          #p .07   r2 .09
summary(lm(t.er.lit~depth,na.action=na.omit,data=data1))          #p .7e-4 r2 .4  
summary(lm(t.er.lit~velocity.mean,na.action=na.omit,data=data1))  #p .99   r2 -.04
summary(lm(t.er.lit~discharge,na.action=na.omit,data=data1))      #p .04   r2 .1
summary(lm(t.er.lit~carbon,na.action=na.omit,data=data1))         #p .99   r2 -.04
summary(lm(t.er.lit~ammonia,na.action=na.omit,data=data1))        #p .6    r2 -.03
summary(lm(t.er.lit~nitrate,na.action=na.omit,data=data1))        #p .4    r2 -.02
summary(lm(t.er.lit~din,na.action=na.omit,data=data1))            #p .5    r2 -.02
summary(lm(t.er.lit~phosphate,na.action=na.omit,data=data1))      #p .5    r2 -.03
summary(lm(t.er.lit~cn.ratio,na.action=na.omit,data=data1))       #p .5    r2 -.02
summary(lm(t.er.lit~canopy,na.action=na.omit,data=data1))         #p .5    r2 -.03
summary(lm(t.er.lit~par.integrative,na.action=na.omit,data=data1))#p .6    r2 -.03
summary(lm(t.er.lit~par.mean,na.action=na.omit,data=data1))       #p .6    r2 -.03
summary(lm(t.er.lit~par.max,na.action=na.omit,data=data1))        #p .9    r2 -.04
summary(lm(t.er.lit~temp.max,na.action=na.omit,data=data1))       #p .3    r2 -.0003
summary(lm(t.er.lit~temp.min,na.action=na.omit,data=data1))       #p .4    r2 -.01 
summary(lm(t.er.lit~temp.mean,na.action=na.omit,data=data1))      #p .4    r2 -.01 
summary(lm(t.er.lit~pr.lit,na.action=na.omit,data=data1))         #p .8    r2 -.04
summary(lm(t.er.lit~cut.capt,na.action=na.omit,data=data1))       #p .5    r2 -.03
summary(lm(t.er.lit~cut.pop,na.action=na.omit,data=data1))        #p .7    r2 -.05
summary(lm(t.er.lit~cut.mass,na.action=na.omit,data=data1))       #p .8    r2 -.05
summary(lm(t.er.lit~cut.mass.m,na.action=na.omit,data=data1))     #p .9    r2 -.06
summary(lm(t.er.lit~sculp.mass.m,na.action=na.omit,data=data1))   #p .008  r2 .3
#pairplot
Z=cbind(data1$slope,data1$stream,data1$depth,data1$pebble,data1$basin,data1$bf,data1$discharge,data1$width,data1$t.er.lit)
colnames(Z)<-c("slope", "stream",    "depth",   "pebble",     "basin",    "bf",     "discharge" ,"width", "t.er.lit")
pairs(Z[,1:9], lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)

lmER1=lm(er.lit~stream+depth+slope+bf,data=data1)
lmER2=lm(er.lit~slope+stream+depth+pebble+basin+bf+discharge+width,data=data1)

summary(lmER1)
lm(formula=lmER1)
anova(lmER1)
drop1(lmER1,test="F")
step(lmER1)

summary(lm(er.lit~basin,na.action=na.omit,data=data1))          #p .02   r2 .2
summary(lm(er.lit~stream,na.action=na.omit,data=data1))         #p .3e-6 r2 .8
summary(lm(er.lit~yr,na.action=na.omit,data=data1))             #p .6    r2 -.03
summary(lm(er.lit~elev,na.action=na.omit,data=data1))           #p .24   r2 .01
summary(lm(er.lit~from.north,na.action=na.omit,data=data1))     #p .9    r2 -.04
summary(lm(er.lit~slope,na.action=na.omit,data=data1))          #p 2e-6  r2 .6 
summary(lm(er.lit~bf,na.action=na.omit,data=data1))             #p .05   r2 .1
summary(lm(er.lit~pebble,na.action=na.omit,data=data1))         #p .05   r2 .1  
summary(lm(er.lit~width,na.action=na.omit,data=data1))          #p .09   r2 .07
summary(lm(er.lit~depth,na.action=na.omit,data=data1))          #p .001  r2 .3  
summary(lm(er.lit~velocity.mean,na.action=na.omit,data=data1))  #p .9    r2 -.04
summary(lm(er.lit~discharge,na.action=na.omit,data=data1))      #p .04   r2 .1
summary(lm(er.lit~carbon,na.action=na.omit,data=data1))         #p .8    r2 -.04
summary(lm(er.lit~ammonia,na.action=na.omit,data=data1))        #p .6    r2 -.03
summary(lm(er.lit~nitrate,na.action=na.omit,data=data1))        #p .6    r2 -.03
summary(lm(er.lit~din,na.action=na.omit,data=data1))            #p .6    r2 -.03
summary(lm(er.lit~phosphate,na.action=na.omit,data=data1))      #p .5    r2 -.03
summary(lm(er.lit~cn.ratio,na.action=na.omit,data=data1))       #p .7    r2 -.04
summary(lm(er.lit~canopy,na.action=na.omit,data=data1))         #p .6    r2 -.03
summary(lm(er.lit~par.integrative,na.action=na.omit,data=data1))#p .4    r2 -.01
summary(lm(er.lit~par.mean,na.action=na.omit,data=data1))       #p .4    r2 -.01
summary(lm(er.lit~par.max,na.action=na.omit,data=data1))        #p .8    r2 -.04
summary(lm(er.lit~temp.max,na.action=na.omit,data=data1))       #p .2    r2 .02
summary(lm(er.lit~temp.min,na.action=na.omit,data=data1))       #p .3    r2 .004 
summary(lm(er.lit~temp.mean,na.action=na.omit,data=data1))      #p .3    r2 .008 
summary(lm(er.lit~pr.lit,na.action=na.omit,data=data1))         #p .8    r2 -.04
summary(lm(er.lit~cut.capt,na.action=na.omit,data=data1))       #p .9    r2 -.05
summary(lm(er.lit~cut.pop,na.action=na.omit,data=data1))        #p .6    r2 -.04
summary(lm(er.lit~cut.mass,na.action=na.omit,data=data1))       #p .7    r2 -.04
summary(lm(er.lit~cut.mass.m,na.action=na.omit,data=data1))     #p .9    r2 -.05
summary(lm(er.lit~sculp.mass.m,na.action=na.omit,data=data1))   #p .004  r2 .4
#pairplot
Z=cbind(data1$stream,data1$slope,data1$depth,data1$basin,data1$bf,data1$pebble,data1$discharge,data1$width,data1$er.lit)
colnames(Z)<-c("stream", "slope",    "depth",   "basin",    "bf",    "pebble",     "discharge" ,"width", "er.lit")
pairs(Z[,1:9], lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)


# predictor transformations
dotchart(data1$slope)
data1$t.slope=log10(1+data1$slope)
dotchart(data1$t.slope)

dotchart(data1$depth)
data1$t.depth=log10(1+data1$depth)
dotchart(data1$t.depth)

dotchart(data1$bf)
data1$t.bf=log10(1+data1$bf)
dotchart(data1$t.bf)
# predictor trans
dotchart(data1$depth)
data1$t.depth=log10(1+data1$depth)
dotchart(data1$t.depth)

dotchart(data1$slope)
data1$t.slope=log10(1+data1$slope)
dotchart(data1$t.slope)

###############################################################################
#using glm to analyze variables
Mer.1=gls(t.er.lit~depth+slope,na.action=na.omit, 
          data=data1, method="ML") #base model
Mer.2=lme(t.er.lit~depth+slope, na.action=na.omit, 
          random = ~1|stream, data=data1, method="ML") #adds random effect
Mer.3=lme(t.er.lit~depth+slope, na.action=na.omit, 
          random = ~1|basin.stream, data=data1, method="ML") #adds nesting

anova(Mer.1, Mer.2, Mer.3)

## # # Analyze residuals
x.er<-data1$t.er.lit[!is.na(data1$t.er.lit)]#removes na values from column
# # Mer.1 residuals
E.1.er<-residuals(Mer.1,type="normalized")
qqnorm(E.1.er)
qqline(E.1.er)
ad.test(E.1.er)
plot(Mer.1) 
plot(x.er, E.1.er)
# # # # M.2 residuals
E.2.er<-residuals(Mer.2,type="normalized")
qqnorm(E.2.er)
qqline(E.2.er)
ad.test(E.2.er)
plot(Mer.2) 
plot(x.er, E.2.er)
abline(0,0)
summary(Mer.2)
# # # M.3 residuals
E.3.er<-residuals(Mer.3,type="normalized")
qqnorm(E.3.er)
qqline(E.3.er)
ad.test(E.3.er)
plot(Mer.3) 
plot(x.er, E.3.er)
#residuals somewhat linear? try alternate variance structures?
#####################################################
#Analyze models with alternate variance structures
#####################################################
#base model from above
Mer.1=gls(t.er.lit~depth+bf, na.action=na.omit, 
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
Mer1.1<-gls(t.er.lit~depth+slope, na.action=na.omit, 
            data=data1, weights=vf1)

Mer1.2<-gls(t.er.lit~depth+slope, na.action=na.omit, 
            data=data1, weights=vf2)

Mer1.3<-gls(t.er.lit~depth+slope, na.action=na.omit, 
            data=data1, weights=vf3)

Mer1.4<-gls(t.er.lit~depth+slope, na.action=na.omit, 
            data=data1, weights=vf4)

Mer1.5<-gls(t.er.lit~depth+slope, na.action=na.omit, 
            data=data1,  weights=varComb(vf1,vf2))

anova(Mer.1,Mer1.1,Mer1.2,Mer1.3,Mer1.5)
anova(Mer.1,Mer1.4)
anova(Mer.1,Mer1.1)

#M1.1
#Analyze alternate variance structure model residuals
E1.1.er<-residuals(Mer1.1,type="normalized")
qqnorm(E1.1.er)
qqline(E1.1.er)
ad.test(E1.1.er)
plot(Mer1.1) 
plot(x.er,E1.1.er)
abline(1,0)
abline(-1,0)
plot(data1$depth,E1.1.er) #have to use na.exclude on model for this to work
plot(data1$bf~E1.1.er)#have to use na.exclude on model for this to work

E1.2.er<-residuals(Mer1.2,type="normalized")
qqnorm(E1.2.er)
qqline(E1.2.er)
ad.test(E1.2.er)
plot(Mer1.2) 
plot(x.er,E1.2.er)

E1.3.er<-residuals(Mer1.3,type="normalized")
qqnorm(E1.3.er)
qqline(E1.3.er)
ad.test(E1.3.er)
plot(Mer1.3) 
plot(x.er,E1.3.er)

E1.4.er<-residuals(Mer1.4,type="normalized")
qqnorm(E1.4.er)
qqline(E1.4.er)
ad.test(E1.4.er)
plot(Mer1.4) 
plot(x.er,E1.4.er)

E1.5.er<-residuals(Mer1.5,type="normalized")
qqnorm(E1.5.er)
qqline(E1.5.er)
ad.test(E1.5.er)
plot(Mer1.5) 
plot(x.er,E1.5.er)

############### ER model validation
op=par(mfrow = c(2, 2)) 
plot(lmER1)
op=par(mfrow = c(2, 2)) 
ER2resid=rstandard(lmER2) 
hist(ER2resid) 
qqnorm(ER2resid) 
plot(ER2resid ~ data1omit$depth, xlab = "depth", ylab = "Residuals",cex.lab=1.5,pch=16) 
abline(0, 0)
plot(ER2resid ~ data1omit$from.north, xlab = "depth", ylab = "Residuals",cex.lab=1.5,pch=16) 
abline(0, 0)

# ER model interpretation
plot(data1omit$depth,data1omit$t.er,cex.lab=1.5,pch=16)
abline(lm(data1omit$t.er~data1omit$depth))
plot(data1omit$from.north,data1omit$t.er,cex.lab=1.5,pch=16)
abline(lm(data1omit$t.er~data1omit$from.north))
