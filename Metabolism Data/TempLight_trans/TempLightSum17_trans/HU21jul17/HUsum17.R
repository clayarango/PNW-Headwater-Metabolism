setwd("N:/Thesis.N/FreshStart_trans/TempLight_trans/TempLightSum17_trans/HU21jul17")
install.packages("chron")
library(chron)
StreamData=read.table("HU21jul17.csv",header=TRUE)  	##time is in seconds since 1970 and in UTC
StreamData$dtime<-chron(StreamData$time/86400)-(7/24)
## Check new chron object dtime and O2 data
## here is where you can target and remove outliers if you want to check before running models 
plot(StreamData$dtime, StreamData$oxy,cex.lab=1.5, cex.axis=1.5, lwd=1,ylim=c(0,9.5), xlab="Time",ylab="Dissolved O2 (mg/L)")
plot(StreamData$dtime, StreamData$oxy,cex.lab=1.5, cex.axis=1.5, lwd=1, xlab="Time",ylab="Dissolved O2 (mg/L)")

# END data important and management
# [2] LOAD O2 SATURATION FUNCTION

osat<- function(temp, bp) {
  
  tstd<-log((298.15-temp) / (273.15 + temp))
  
  a0<-2.00907
  a1<-3.22014
  a2<-4.0501
  a3<-4.94457
  a4<- -0.256847
  a5<- 3.88767
  
  u<-10^(8.10765-(1750.286/(235+temp)))
  
  sato<-(exp(a0 + a1*tstd + a2*tstd^2 + a3*tstd^3 + a4*tstd^4+a5*tstd^5))*((bp-u)/(760-u))*1.42905
  sato
}

#end of function
# END loading O2 SAT calc #
# [3] LOAD BAROMETRIC PRESSURE FUNCTION
## Function to correct barometric pressure for altitude. From Colt (2012).
## This function gives bp in mmHg for altitude given nearby measurement of standardized barometric pressure. 
## temp is degC 
## alt is m
## bpst is in inches of Hg and the sort of data you will find from U.S. weather stations.  Delete the *25.4 if you in the rest of the world where pressure is given in metric units
####function returns mm of Hg
bpcalc<- function(bpst, alt) {
  bpst*25.4*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+15)))
}
#end of function
# END loading BP calc function ####### 
# [4] LOAD GAS EXCHANGE (K) FUNCTIONS
# NOTE: The functions you use will depend on the methods used to estimate K (O2, propane, SF6). The temperature correction (Kcor) is embedded in the models below, but the Kcor function must be loaded into R before running the model.
# UNITS are day^(-1)
## This code does the opposite of Kcor below; it estimates K600 for KO2 at a given temperature. From Wanninkhof (1992).
K600fromO2<-function (temp, KO2) {
  ((600/(1800.6 - (120.1 * temp) + (3.7818 * temp^2) - (0.047608 * temp^3)))^-0.5) * KO2
}
#end of function
# This calculates K600 from K measured in a propane addition. From Raymond et al. (2012).
K600frompropane<-function (temp, Kpropane) {
  ((600/(2864 - (154.14 * temp) + (3.791 * temp^2) - (0.0379 * temp^3)))^-0.5) * Kpropane
}
#end of function
# This calculates K600 from K measured in a SF6 addition. From Raymond et al. (2012).
K600fromSF6<-function(temp, KSF6) {
  ((600/(3255.3-(217.13*temp)+(6.837*temp^2)-(0.08607*temp^3)))^-0.5)*KSF6
}
#end of function
# This function calculates KO2 at any given tempertaure from K600. via schmidt number scaling.  The scaling equation if From JÃ¤hne et al. (1987), Schmidt number conversions from Wanninkhof et al. 1992.
Kcor<-function (temp,K600){
  K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}
#end of function
# END loading K functions
# [5] NIGHTTIME REGRESSION to estimate K -- OPTIONAL
### nighttime regression code to estimate K. Approach as in Hornberger and Kelly (1975).
## o2 file is your oxygen data (defined in subsetting)
## bp is barometric pressure in mm Hg for your site, 
## ts is the time step in MINUTES (not days as in metabolism code below)

nightreg<-function(o2file, bp, ts){
  
  temp<-o2file$temp
  oxy<-o2file$oxy
  
  # moving average on oxy data
  oxyf1<-filter(o2file$oxy,rep(1/3,3), sides=2)
  
  # trim the ends of the oxy data
  oxyf2<- oxyf1[c(-1,-length(oxyf1))]
  
  # calculate delO/delt; convert to units of days by ts in min-1*1440
  deltaO2<-((oxyf2[-1]-oxyf2[-length(oxyf2)])/ts)*1440
  
  # Trim the first two and last one from the temp data to match the filter oxy data
  temptrim<-temp[c(-2:-1,-length(temp))]
  
  # calc the dodef
  satdef<-osat(temptrim,bp)-oxyf2[-1]
  
  # fit linear model and plot using linear model fitting (lm) and abline functions in R stats package
  nreg<-lm(deltaO2~satdef)
  plot(satdef,deltaO2)
  abline(nreg)
  
  # use coef function in R stats package to get lm coefficients
  coeff<-coef(nreg)
  
  # output gives lm coeff and K600 (converted from nighttime regression estimate of KO2)
  out<-list(coeff, K600fromO2(mean(temp), coeff[2]))
  out
  
}
#end of function
# NOTE: this approach works better for some sites/dates than others; always check that your model fit is good and that your K600 estimate makes sense!
#Call as:  The first argument in the function defines when to pull data.  In this case on 10/27/204 (for spring creek) between 18:05 and 23:00

nightreg(StreamData[StreamData$dtime>=as.numeric(chron(dates="07/19/17", times="18:10:00")) & StreamData$dtime<=as.numeric(chron(dates="07/21/17", times="10:00:00")), ], bp=660, ts=10)
# END nighttime regression calculation of K #
# check that light makes sense
plot(StreamData$dtime, StreamData$light,cex.lab=1.5, cex.axis=1.5, lwd=1, xlab="Time",ylab="PAR (umol photons/s/m2)")
# END light
# [7] MODEL v1 - RIVERMETABK; solves for GPP, ER, AND K
##
# [7a] LOAD RIVERMETABK functions
###########################################

## VERSION 1 of possible metabolism model

##### now estimate metabolism

## parameter names (and units) for BOTH rivermetab functions (solving for or fixing K600):
## oxy.mod = modeled O2 (mg/L)
## oxy = O2 data (mg/L)
## MET[1] = GPP = estimated daily gross primary production (g O2 m-2 d-1); + flux of O2 production
## MET[2] = ER = estimated daily ecosystem respiration (g O2 m-2 d-1); - flux of O2 consumption 
## MET[3] = K = estimated K600 (d-1)
## z = estimated site depth (m)
## Kcor = KO2 corrected for temperature, calculated from K600 (d-1)
## bp = barometric pressure (mmHg)
## temp = water temperature (C)
## light = PAR data (or modlight from above light function)
## ts = time step of O2 data from logger (10 min intervals --> units are day, so 10/1440=0.06944)

# Model to calculate GPP, ER and K simultaneously
# This model is advantageous in the sense that one can estimate K from the data, but beware that it may be an overparameterized model.  
# Note that you have the opportunity below to use your light data or modeled light (here as data$light estimated for spring creek and french creek with the light model function above). You decide and modify the function and data names accordingly.

# rivermetabK is the master function that calls the MLE function (onestationmleK) and plotting function (onestationplot) below
rivermetabK<-function(o2file, z, bp, ts){
  
  ##pull data out of loaded o2file (subset in CALL function) and give it a name. 
  temp<-o2file$temp
  oxy<-o2file$oxy
  light<-o2file$light
  
  ##This calls onestationmleK (below) to calculate metabolism by non-linear minimization. We use nlm() function in R stats package; for more information see "help(nlm)"  The first argument is the function to be minimized, the second defines the starting values.  The function that is minimized (onestationmleK, see below) always has as its first argument, a vector of parameters (MET), and second argument a vector of starting values for those parameters, p=c(3,-5, 10).
  
  river.mle<-nlm(onestationmleK, p=c(3,-5, 10), oxy=oxy, z=z,temp=temp,light=light, bp=bp, ts=ts)
  
  ##plot modeled and measaured O2 given MLE estimates of GPP, ER, and K600.  It calls a function below onestationplot()
  
  onestationplot(GPP=river.mle$estimate[1], ER=river.mle$estimate[2],oxy=oxy,z=z,temp=temp,light=light, K=river.mle$estimate[3], bp=bp, ts=ts)
  
  ##return GPP, ER, K600, and MLE values
  b <- list(GPP=river.mle$estimate[1], ER=river.mle$estimate[2], K600=river.mle$estimate[3], neglogL=river.mle$minimum[1])
  b
}
# end of function
# This function returns the negative log likelihood value given O2 data and estimates of GPP, ER, and K (which is vector MET); is included in master rivermetabK function above
onestationmleK<-function(MET,temp, oxy, light, z, bp, ts) {
  
  # create new vector for modeled O2
  oxy.mod<-numeric(length(data))
  # give starting value from oxygen data; this is the only time O2 data is used to model GPP and ER
  oxy.mod[1]<-oxy[1]
  
  # this is the metabolism equation as in Van de Bogert et al 2007 L&OMethods
  for (i in 2:length(oxy)) {oxy.mod[i]<-oxy.mod[i-1]+((MET[1]/z)*(light[i]/sum(light)))+ MET[2]*ts/z+(Kcor(temp[i],MET[3]))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }
  
  ##below is MLE calculation; output is -log likelihood
  # diff in measured and modeled O2
  sqdiff<-(oxy-oxy.mod)^2 
  # likelihood function
  L <- length(oxy) * (log(((sum(sqdiff)/length(oxy))^0.5))+0.5*log(6.28)) + ((2*sum(sqdiff)/length(oxy))^-1) * sum(sqdiff)
  L
}
# end of function
# this function plots modeled O2 and O2 data from estimates of daily GPP, ER, and K; is included in master rivermetabK function above
# Calls same metabolism equation as in mle function, but plots modeled O2 as a function of GPP, ER, and K estimates from mle
# use this to visually assess your model estimates (should be good agreement between modeled and measured O2)
onestationplot<-function(GPP, ER, oxy, z, temp, K, light, bp, ts) {
  
  oxy.mod<-numeric(length(oxy))
  oxy.mod[1]<-oxy[1]
  
  # this is the metabolism equation as in Van de Bogert et al (2007) L&OMethods
  for (i in 2:length(oxy)) { oxy.mod[i]<-oxy.mod[i-1]+((GPP/z)*(light[i]/sum(light)))+ ER*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }
  
  plot(seq(1:length(oxy)),oxy.mod, type="l",xlab="Time", ylab="Dissolved oxygen  (mg/L)", cex.lab=1.5, cex.axis=1.5, lwd=2 )
  points(seq(1:length(oxy)),oxy)
}
# end of function

####### END LOADING rivermetabK function #
# [7b] CALL RIVERMETABK function
###########################################

## Call as:  z is depth in m, bp is im mmHg, ts is time steps in days.

# for spring creek data; 10/28/14
###########################################

## Call as:  z is depth in m, bp is im mmHg, ts is time steps in days.

# for spring creek data; 10/28/14
rivermetab<-function(o2file, z, bp, ts, K){
  
  ##pull data out of loaded o2file (subset in CALL function) and give it a name. 
  temp<-o2file$temp
  oxy<-o2file$oxy
  light<-o2file$light
  
  ##calculate metabolism by non linear minimization of MLE function (below)
  river.mle<-nlm(onestationmle, p=c(3,-5), oxy=oxy, z=z,temp=temp,light=light, bp=bp, ts=ts, K=K)
  
  ##plot data; uses same plot function as given for rivermetabK above (relisted below to keep each model as own unit)
  onestationplot(GPP=river.mle$estimate[1], ER=river.mle$estimate[2],oxy=oxy,z=z,temp=temp,light=light, K=K, bp=bp, ts=ts)
  
  ##return GPP, ER, and MLE value
  b<-list(GPP= river.mle$estimate[1], ER= river.mle$estimate[2],  neglogL= river.mle$minimum[1])
  b
}
# end of function

# function returns the likelihood value given O2 data and estimates of GPP, ER (which is vector MET); is included in master rivermetab function above; K600 is fixed
onestationmle<-function(MET,temp, oxy, light, z, bp, ts, K) {
  
  oxy.mod<-numeric(length(data))
  oxy.mod[1]<-oxy[1]
  
  # this is the metabolism equation as in Van de Bogert et al 2007 L&OMethods
  for (i in 2:length(oxy)) {oxy.mod[i]<-oxy.mod[i-1]+((MET[1]/z)*(light[i]/sum(light)))+ MET[2]*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }
  
  ## below is MLE calculation; output is likelihood
  # diff in measured and modeled O2
  sqdiff<-(oxy-oxy.mod)^2 
  # likelihood function
  L <- length(oxy)*(log(((sum(sqdiff)/length(oxy))^0.5)) +0.5*log(6.28)) + ((2*sum(sqdiff)/length(oxy))^-1)*sum(sqdiff)
  L
}
#end of function

# (As in rivermetabK)
# this function plots modeled O2 and O2 data from estimates of daily GPP and ER; is included in master rivermetab function above
# Calls same metabolism equation as in mle function, but plots modeled O2 as a function of GPP, ER estimates from mle
# use this to visually assess your model estimates (should be good agreement between modeled and measured O2)
onestationplot<-function(GPP, ER, oxy, z, temp, K, light, bp, ts) {
  
  oxy.mod<-numeric(length(oxy))
  oxy.mod[1]<-oxy[1]
  
  # this is the metabolism equation as in Van de Bogert et al (2007) L&OMethods
  for (i in 2:length(oxy)) { oxy.mod[i]<-oxy.mod[i-1]+((GPP/z)*(light[i]/sum(light)))+ ER*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }
  
  plot(seq(1:length(oxy)),oxy.mod, type="l",xlab="Time", ylab="Dissolved oxygen  (mg/L)", cex.lab=1.5, cex.axis=1.5, lwd=2 )
  points(seq(1:length(oxy)),oxy)
  
  print(summary(lm(oxy~oxy.mod)))
  
}
# end of function


####### END loading rivermetab function ####### 
## Call as:  z is depth in m, bp is im mmHg, ts is time steps in days.

# 10 MIN timestep!! literature k600 values based on slope
rivermetab(o2file=StreamData[StreamData$dtime>=as.numeric(chron(
  dates="07/21/17", times="18:22:00")) & StreamData $dtime<=as.numeric(chron(
    dates="07/23/17", times="05:32:00")), ], 
  z=0.046, 
  bp=724, 
  ts=0.006944, 
  K=80.768051)
# 10 MIN timestep!! my derived k600 values based on velocity
rivermetab(o2file=StreamData[StreamData$dtime>=as.numeric(chron(
  dates="07/21/17", times="18:22:00")) & StreamData $dtime<=as.numeric(chron(
  dates="07/23/17", times="05:32:00")), ], 
  z=0.046, 
  bp=724, 
  ts=0.006944, 
  K=20.0675500)
# 10 MIN timestep!! full model as is 
rivermetabK(o2file=StreamData[ StreamData $dtime>=as.numeric(chron(
  dates="07/21/17", times="18:22:00")) & StreamData $dtime<=as.numeric(chron(
  dates="07/23/17", times="05:32:00")), ], 
  z=0.046, 
  bp=724, 
  ts=0.006944)

