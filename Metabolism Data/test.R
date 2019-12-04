### Master script to analyze PNW Metabolism data via streammetabolizer
### Created by AJR
### Created on 2019-12-04


### Read in data from git repo
test.data<-read.csv(file="./Metabolism Data/test.csv")

#Install stream metabolizer:
install.packages("streamMetabolizer", dependencies=TRUE, 
                repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))

#read in packages
library(chron)
library(streametabolizer)

# LOAD O2 SATURATION FUNCTION

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

test.data$DO.sat<-osat(test.data$temp.water, test.data$bp)

