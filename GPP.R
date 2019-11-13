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


##################### model selections

# GPP transformation
dotchart(data1$gpp.lit)
ad.test(data1$gpp.lit)
data1$t.gpp.lit=log(1+data1$gpp.lit)^(1/2)
ad.test(data1$t.gpp.lit)
hist(data1$t.gpp.lit,xlab="transfromed gpp",cex.lab=1.5)
dotchart(data1$t.gpp.lit,xlab="transformed gpp",cex.lab=1.5,pch=16)
# GPP LM models
summary(lm(gpp.lit~season,na.action=na.omit,data=data1))         #p .3    r2 .01 
summary(lm(gpp.lit~season.yr,na.action=na.omit,data=data1))      #p .12   r2 .09 
summary(lm(gpp.lit~basin,na.action=na.omit,data=data1))          #p .08   r2 .13 
summary(lm(gpp.lit~stream,na.action=na.omit,data=data1))         #p .05   r2 .4  
summary(lm(gpp.lit~yr,na.action=na.omit,data=data1))             #p .3    r2 .004
summary(lm(gpp.lit~elev,na.action=na.omit,data=data1))           #p .09   r2 .08  
summary(lm(gpp.lit~from.north,na.action=na.omit,data=data1))     #p .8    r2 -.04
summary(lm(gpp.lit~slope,na.action=na.omit,data=data1))          #p .08   r2 .09 
summary(lm(gpp.lit~bf,na.action=na.omit,data=data1))             #p .9    r2 -.04
summary(lm(gpp.lit~pebble,na.action=na.omit,data=data1))         #p .0001 r2 .45
summary(lm(gpp.lit~width,na.action=na.omit,data=data1))          #p .15   r2 .05
summary(lm(gpp.lit~depth,na.action=na.omit,data=data1))          #p .04   r2 .1
summary(lm(gpp.lit~velocity.mean,na.action=na.omit,data=data1))  #p .2    r2 .03
summary(lm(gpp.lit~discharge,na.action=na.omit,data=data1))      #p .9    r2 -.04
summary(lm(gpp.lit~carbon,na.action=na.omit,data=data1))         #p .9    r2 -.04
summary(lm(gpp.lit~ammonia,na.action=na.omit,data=data1))        #p .1    r2 .06
summary(lm(gpp.lit~nitrate,na.action=na.omit,data=data1))        #p .4    r2 -.01
summary(lm(gpp.lit~din,na.action=na.omit,data=data1))            #p .4    r2 -.005
summary(lm(gpp.lit~phosphate,na.action=na.omit,data=data1))      #p .16   r2  .04
summary(lm(gpp.lit~cn.ratio,na.action=na.omit,data=data1))       #p .5    r2  -.02
summary(lm(gpp.lit~canopy,na.action=na.omit,data=data1))         #p .1    r2  .07
summary(lm(gpp.lit~par.integrative,na.action=na.omit,data=data1))#p .26   r2  .01
summary(lm(gpp.lit~par.mean,na.action=na.omit,data=data1))       #p .26   r2  .01
summary(lm(gpp.lit~par.max,na.action=na.omit,data=data1))        #p .8    r2  -.04
summary(lm(gpp.lit~temp.max,na.action=na.omit,data=data1))       #p .07   r2  .1  
summary(lm(gpp.lit~temp.min,na.action=na.omit,data=data1))       #p .07   r2  .09 
summary(lm(gpp.lit~temp.mean,na.action=na.omit,data=data1))      #p .08   r2  .08 
summary(lm(gpp.lit~cut.capt,na.action=na.omit,data=data1))       #p .7    r2 -.05
summary(lm(gpp.lit~cut.pop,na.action=na.omit,data=data1))        #p .5    r2 -.03
summary(lm(gpp.lit~cut.mass,na.action=na.omit,data=data1))       #p .2    r2 .04
summary(lm(gpp.lit~cut.mass.m,na.action=na.omit,data=data1))     #p .6    r2 -.03
summary(lm(gpp.lit~sculp.mass.m,na.action=na.omit,data=data1))   #p .02   r2 .2
summary(lmGPP)
#pairplot non transformed
Z=cbind(data1$basin,data1$stream,data1$elev,data1$slope,data1$pebble,data1$depth,data1$temp.max,data1$gpp.lit)
colnames(Z)<-c("basin", "stream",    "elev",    "slope",    "pebble",   "depth",   "temp.max",    "gpp.lit")
pairs(Z[,1:8], lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)


# do transformed ones!!
summary(lm(t.gpp.lit~season,na.action=na.omit,data=data1))         #p .2   r2 .02
summary(lm(t.gpp.lit~season.yr,na.action=na.omit,data=data1))      #p .1   r2 .1 
summary(lm(t.gpp.lit~basin,na.action=na.omit,data=data1))          #p .1   r2 .1 
summary(lm(t.gpp.lit~stream,na.action=na.omit,data=data1))         #p .03  r2 .4 
summary(lm(t.gpp.lit~yr,na.action=na.omit,data=data1))             #p .3   r2 -.002
summary(lm(t.gpp.lit~elev,na.action=na.omit,data=data1))           #p .06  r2 .1  
summary(lm(t.gpp.lit~from.north,na.action=na.omit,data=data1))     #p .9   r2 -.04
summary(lm(t.gpp.lit~slope,na.action=na.omit,data=data1))          #p .06  r2 .1 
summary(lm(t.gpp.lit~bf,na.action=na.omit,data=data1))             #p .7   r2 -.03.
summary(lm(t.gpp.lit~pebble,na.action=na.omit,data=data1))         #p .0001 r2 .4
summary(lm(t.gpp.lit~width,na.action=na.omit,data=data1))          #p .2   r2 .02
summary(lm(t.gpp.lit~depth,na.action=na.omit,data=data1))          #p .04  r2 .1
summary(lm(t.gpp.lit~velocity.mean,na.action=na.omit,data=data1))  #p .1   r2 .07
summary(lm(t.gpp.lit~discharge,na.action=na.omit,data=data1))      #p .7   r2 -.04
summary(lm(t.gpp.lit~carbon,na.action=na.omit,data=data1))         #p .5   r2 -.02
summary(lm(t.gpp.lit~ammonia,na.action=na.omit,data=data1))        #p .09  r2 .08
summary(lm(t.gpp.lit~nitrate,na.action=na.omit,data=data1))        #p .5   r2 -.02
summary(lm(t.gpp.lit~din,na.action=na.omit,data=data1))            #p .4   r2 -.02
summary(lm(t.gpp.lit~phosphate,na.action=na.omit,data=data1))      #p .2   r2  .03
summary(lm(t.gpp.lit~cn.ratio,na.action=na.omit,data=data1))       #p .2   r2  .02
summary(lm(t.gpp.lit~canopy,na.action=na.omit,data=data1))         #p .1   r2  .05
summary(lm(t.gpp.lit~par.integrative,na.action=na.omit,data=data1))#p .3   r2  .02
summary(lm(t.gpp.lit~par.mean,na.action=na.omit,data=data1))       #p .3   r2  .01
summary(lm(t.gpp.lit~par.max,na.action=na.omit,data=data1))        #p .9   r2  -.04
summary(lm(t.gpp.lit~temp.max,na.action=na.omit,data=data1))       #p .06  r2  . 1
summary(lm(t.gpp.lit~temp.min,na.action=na.omit,data=data1))       #p .06  r2  .1 
summary(lm(t.gpp.lit~temp.mean,na.action=na.omit,data=data1))      #p .07  r2  .09
summary(lm(t.gpp.lit~cut.capt,na.action=na.omit,data=data1))       #p .5   r2 -.03
summary(lm(t.gpp.lit~cut.pop,na.action=na.omit,data=data1))        #p .7   r2 -.04
summary(lm(t.gpp.lit~cut.mass,na.action=na.omit,data=data1))       #p .2   r2 .04
summary(lm(t.gpp.lit~cut.mass.m,na.action=na.omit,data=data1))     #p .7   r2 -.04
summary(lm(t.gpp.lit~sculp.mass.m,na.action=na.omit,data=data1))   #p .02  r2 .2

#pairplot transformed  Change to appropriate variables
Z=cbind(data1$pebble,data1$stream,data1$depth,data1$elev,data1$slope,data1$temp.max,data1$season.yr,data1$basin,data1$t.gpp.lit)
colnames(Z)<-c("pebble", "stream",    "depth",    "elev",    "slope",   "temp.max",   "season.yr",   "basin", "t.gpp.lit")
pairs(Z[,1:9], lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)




ad.test(residuals(lmGPP))
hist(rstandard(lmGPP)) 
qqnorm(rstandard(lmGPP)) 
#    

lmGPP1=lm(t.gpp.lit~stream+season.yr+par.integrative+din,data=data1)
summary(lmGPP1)

lmGPP1.2=lm(t.gpp.lit~pebble+stream+season.yr,data=data1)
summary(lmGPP1.2)

lmGPP1.3=lm(t.gpp.lit~stream+season.yr,data=data1)
summary(lmGPP1.3)

lmGPP2=lm(t.gpp.lit~pebble+stream+depth+elev+slope+temp.max+season.yr+basin,data=data1)
summary(lmGPP1)
lm(formula=lmGPP1)
anova(lmGPP1)
drop1(lmGPP1,test="F")
step(lmGPP1)
plot(lmGPP1)
plot(rstandard(lmGPP1)~(data1$t.gpp.lit[!is.na(data1$t.gpp.lit)]))
abline(1,0)
abline(-1,0)

#transform explanatory
dotchart(data1$pebble)
data1$t.pebble=log10(1+data1$pebble)
dotchart(data1$t.pebble)

dotchart(data1$depth)
data1$t.depth=log10(1+data1$depth)
dotchart(data1$t.depth)

dotchart(data1$temp.max)
data1$t.temp.max=(1+data1$temp.max)^2
dotchart(data1$t.temp.max)
ad.test(data1$t.temp.max)
###############################################################################
#using glm to analyze variables
M.1=gls(t.gpp.lit~season.yr+depth,na.action=na.omit, 
        data=data1, method="ML") #base model
M.2=lme(t.gpp.lit~season.yr+depth, na.action=na.omit, 
        random = ~1|stream, data=data1, method="ML") #adds random effect

M.3=lme(t.gpp.lit~season.yr+canopy, na.action=na.omit, 
        random = ~1|basin.stream, data=data1, method="ML") #adds nesting

anova(M.1, M.2, M.3)

#Analyze residuals
x.gpp<-data1$t.gpp.lit[!is.na(data1$t.gpp.lit)]#removes na values from column
x.pebble=data1$t.pebble[!is.na(data1$gpp.lit)]
x.depth=data1$t.depth[!is.na(data1$gpp.lit)]
x.temp.max=data1$t.temp.max[!is.na(data1$gpp.lit)]
# M.1 residuals
E.1.gpp<-residuals(M.1,type="normalized")
qqnorm(E.1.gpp)
qqline(E.1.gpp)
ad.test(E.1.gpp)
plot(M.1) 
plot(x.gpp, E.1.gpp)
summary(lm(x.gpp~E.1.gpp))
plot(data1$gpp.lit,data1$par.max)
# M.2 residuals
E.2.gpp<-residuals(M.2,type="normalized")
qqnorm(E.2.gpp)
qqline(E.2.gpp)
ad.test(E.2.gpp)
plot(M.2) 
plot(x.gpp, E.2.gpp)
summary(lm(x.gpp~E.2.gpp))

acf(E.2.gpp, na.action=na.pass,
    main="Auto-correlation plot for residuals")#check for autocorrelation in residuals

# M.3 residuals
E.3.gpp<-residuals(M.3,type="normalized")
qqnorm(E.3.gpp)
qqline(E.3.gpp)
ad.test(E.3.gpp)
plot(M.3) 
plot(x.gpp, E.3.gpp)
summary(lm(x.gpp~E.3.gpp))



#residuals somewhat linear? try alternate variance structures?

#####################################################
#Analyze models with alternate variance structures
#####################################################

#base model from above
M.2=lme(t.gpp.lit~depth+temp.max, na.action=na.omit, 
        random = ~1|stream, data=data1) #run with method="REML" default for comparison

#alternate variance structures
vf1=varIdent(form = ~1|stream)
vf2=varPower(form = ~fitted(.))
vf3=varExp(form = ~fitted(.))
vf4=varConstPower(form = ~fitted(.))
#vf5=varPower(form = ~fitted (.)|site)
#vf6=varExp(form = ~fitted(.)|site)
#vf7=varConstPower(form = ~fitted(.)|site)
#vf8=varIdent(form = ~1|f.sample.event)

#alternate models+temp.max
M1.1<-lme(t.gpp.lit~season.yr+depth, na.action=na.omit, 
          random = ~1|stream, data=data1, weights=vf1)

M1.2<-lme(t.gpp.lit~depth+temp.max, na.action=na.omit, 
          random = ~1|stream, data=data1, weights=vf2)

M1.3<-lme(t.gpp.lit~depth+temp.max, na.action=na.omit, 
          random = ~1|stream, data=data1, weights=vf3)

M1.4<-lme(t.gpp.lit~depth+temp.max, na.action=na.omit, 
          random = ~1|stream, data=data1, weights=vf4)

M1.5<-lme(t.gpp.lit~depth+temp.max, na.action=na.omit, 
          random = ~1|stream, data=data1,  weights=varComb(vf1,vf2))

anova(M.2,M1.1,M1.2,M1.3,M1.4,M1.5)
anova(M.1,M1.1)
#M1.3, expoential variance of fit

#Analyze alternate variance structure model residuals

E1.1.gpp<-residuals(M1.1,type="normalized")
qqnorm(E1.1.gpp)
qqline(E1.1.gpp)
ad.test(E1.1.gpp)
plot(M1.1) 
plot(x.gpp, E1.1.gpp)
summary(lm(x.gpp~E1.1.gpp))

plot(x.pebble, E1.1.gpp)
plot(x.depth, E1.1.gpp)
plot(x.temp.max, E1.1.gpp)

E1.2.gpp<-residuals(M1.2,type="normalized")
qqnorm(E1.2.gpp)
qqline(E1.2.gpp)
ad.test(E1.2.gpp)
plot(M1.2) 
plot(x.gpp, E1.2.gpp)
summary(lm(x.gpp~E1.2.gpp))

E1.3.gpp<-residuals(M1.3,type="normalized")
qqnorm(E1.3.gpp)
qqline(E1.3.gpp)
ad.test(E1.3.gpp)
plot(M1.3) 
plot(x.gpp, E1.3.gpp)
summary(lm(x.gpp~E1.3.gpp))


E1.4.gpp<-residuals(M1.4,type="normalized")
qqnorm(E1.4.gpp)
qqline(E1.4.gpp)
ad.test(E1.4.gpp)
plot(M1.4) 
plot(x.gpp, E1.4.gpp)
summary(lm(x.gpp~E1.4.gpp))


E1.5.gpp<-residuals(M1.5,type="normalized")
qqnorm(E1.5.gpp)
qqline(E1.5.gpp)
ad.test(E1.5.gpp)
plot(M1.5) 
plot(x.gpp, E1.5.gpp)
summary(lm(x.gpp~E1.5.gpp))

##################################################################################
# GPP model validation
op=par(mfrow = c(2, 2)) 
plot(lmGPP1)
op=par(mfrow = c(2, 2)) 
GPP1resid=rstandard(lmGPP1) 

hist(GPP2resid) 
qqnorm(GPP2resid) 
plot(GPP2resid ~ data1omit$pebble, xlab = "pebble", ylab = "Residuals",cex.lab=1.5,pch=16) 
abline(0, 0)
plot(GPP2resid ~ data1omit$from.north, xlab = "aspect from N", ylab = "Residuals",cex.lab=1.5,pch=16) 
abline(0, 0)

#   plot(x,y) but abline(lm(y~x)) !!!!
# GPP model interpretation
plot(data1$pebble,data1$t.gpp.lit,cex.lab=1.5,pch=16)
abline(lm(data1$t.gpp~data1$pebble))
summary(lm(data1$t.gpp~data1$pebble))
legend("topright", legend=c("R2= 0.45", "P= 1e-4"),bty="n")

plot(data1$temp.max,data1$t.gpp.lit,cex.lab=1.5,pch=16)
abline(lm(data1$t.gpp.lit~data1$temp.max))
summary(lm(data1$t.gpp.lit~data1$temp.max))
legend("topright", legend=c("R2= 0.10", "P=0.07"),bty="n")



