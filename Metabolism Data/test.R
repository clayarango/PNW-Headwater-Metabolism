### Master script to analyze PNW Metabolism data via streammetabolizer
### Created by AJR
### Created on 2019-12-04

### AJR got the base model and a model with a constrained K600 (based on slope) prior running
### AJR next steps: work on the calc_light compared to measured light models, but calc_light 
### doesn't seem to give the right data (i.e., light peaks at ~10PM based on calc_light right now...)


### Read in data from git repo
test.data<-read.csv(file="./Metabolism Data/test.csv")

#If you don't have streamMetabolizer installed, first install it by uncommenting the next two lines
#install.packages("streamMetabolizer", dependencies=TRUE, 
#                repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))

#If you already have streamMetabolizer installed, make sure you are using the most up-to-date version by uncommenting the next two lines
update.packages(oldPkgs=c("streamMetabolizer", "unitted"), depedencies=T, 
                repos=c("https://owi.usgs.gov/R", "https://cran.rstudio.com"))


#read in packages
library(chron)
library(streamMetabolizer)
library(dplyr)
library(tidyr)
library(ggplot2)


### Convert time
test.data$datetime<-as.POSIXct(test.data$time, origin="1970-01-01", tz='America/Los_Angeles') # time is in seconds since 1970, this converts to %m/%D/%Y %HH:%MM:%SS in the local tz of the sites
#lubridate::tz(test.data$datetime) # Check to make sure you get it into the appropriate tz
#Note, solar time calculations need to account for each individual site's longitude - need to check if this can reference another column or not
test.data$solar.time<-calc_solar_time(test.data$datetime, longitude=120.782771) # calculate solar time which is needed for streammetabolizer

# Calculate O2 saturation based on temperature and bp
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
#Now run it on the data to calulcate DO.sat
test.data$DO.sat<-osat(test.data$temp.water, test.data$bp)


#Check that DO %sat seemingly makes sense:
plot(test.data$datetime, (test.data$DO.obs/test.data$DO.sat)*100)

# Try having stream metabolizer get light for you instead of using canopy-influenced light
test.data$calc.light<-calc_light(test.data$solar.time, latitude=47.296801, longitude=120.705743)

#Compare measured PAR with streammetabolizer PAR:
plot(test.data$datetime, test.data$calc.light)
plot(test.data$datetime, test.data$light)
plot(test.data$datetime, test.data$solar.time)


#Look at the data through streammetabolizer's eyes:
windows(width=7, height=3)
test.data %>% unitted::v() %>%
  mutate(DO.pctsat=100*(DO.obs/DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() + scale_color_discrete('variable')

data<-test.data[,9:14]
light.test.data<-test.data[,9:14]
light.test.data$light<-test.data$calc.light

#Pick a model in stream metabolizer:
base.model<-mm_name(
  type="bayes", pool_K600='none')

base.model_specs<-specs(base.model, day_start=10, day_end=34, burnin_steps=100, saved_steps=100)
k.slope.model_specs<-specs(base.model, day_start=10, day_end=34, K600_daily_meanlog=log(80.768051), K600_daily_sdlog=0.01, burnin_steps=100, saved_steps=100)


base.fit<-metab(base.model_specs, data=data)
base.light.fit<-metab(base.model_specs, data=light.test.data)
k.slope.fit<-metab(k.slope.model_specs, data=data)
k.slope.light.fit<-metab(k.slope.model_specs, data=light.test.data)

predict_metab(base.fit)
predict_metab(k.slope.fit)
predict_metab(k.slope.light.fit)
plot_metab_preds(k.slope.fit)
plot_metab_preds(k.slope.light.fit)
plot_DO_preds(k.slope.fit)
plot_DO_preds(k.slope.light.fit)
plot_DO_preds(base.light.fit)

predict_DO(k.slope.fit)
base.mcmc<-get_mcmc(base.fit)
slope.mcmc<-get_mcmc(k.slope.fit)
traceplot(base.mcmc, pars="K600_daily")
traceplot(slope.mcmc, pars="ER")
get_fit(base.fit)$overall %>%
  select(ends_with('Rhat'))
