install.packages("tidyverse")
install.packages("nlme")
install.packages("nortest")
install.packages("multcomp")
install.packages("MASS")
install.packages("gridExtra")
install.packages("vegan")
install.packages("agricolae")

library(tidyverse)
library(nlme)
library(nortest)
library(multcomp)
library(MASS)
library(gridExtra)
library(vegan)
library(agricolae)

data1=read.csv("all_trans_short_r.csv",header=TRUE) #load data with NA's
data1cut=read.csv("trans_short_cut_r.csv",header=TRUE)
dataint=read.csv("interaction2by3.csv",header=TRUE)


data1$t.gpp.lit=log(1+data1$gpp.lit)^(1/2)
data1$t.er.lit=log(1+abs(data1$er.lit))
data1$t.cut.mass.m=log(1+data1$cut.mass.m)

# y axis title help =      (bquote('PAR ('*mu~ 'mol' ~O[2]~ m^-2~d^-1*')'))
#  annotation help (bottom of page)  https://rstudio-pubs-static.s3.amazonaws.com/136237_170402e5f0b54561bf7605bdea98267a.html
    #  annotate("text", x=0.04, y=.75,  size=3, hjust=0,parse = TRUE, label = as.character(expression(paste(R^"2","=0.13"))))

#############  Carbon
xcarbon=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(carbon.mean = mean(carbon, na.rm = TRUE), # na.rm = TRUE to remove missing values
            carbon.sd=sd(carbon, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            carbon.n = sum(!is.na(carbon)), # of observations, excluding NAs. 
            carbon.se=carbon.sd/sqrt(carbon.n))

plot(HSD.test(aov(data1$carbon~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Carbon")
anova(aov(data1$carbon~data1$stream+data1$season.yr))

ggplot(data=xcarbon,aes(x=season.yr,y=carbon.mean))+geom_bar(stat="identity")

ggplot(data=xcarbon, 
       aes(x=season.yr, y=carbon.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=carbon.mean, ymax=carbon.mean+carbon.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Carbon ('*'DOC'~~ mg~L^-1*')'))+
  ylim(0,10) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=4.6, label="B", size=3, hjust=0) +
  annotate("text", x=1.95, y=5.4, label="B", size=3, hjust=0)+
  annotate("text", x=2.95, y=9.2, label="A", size=3, hjust=0)+
  annotate("text", x=2.8, y=10, label="P=0.00012", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.carbon1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############  Phosphate
xphosphate=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(phosphate.mean = mean(phosphate, na.rm = TRUE), # na.rm = TRUE to remove missing values
            phosphate.sd=sd(phosphate, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            phosphate.n = sum(!is.na(phosphate)), # of observations, excluding NAs. 
            phosphate.se=phosphate.sd/sqrt(phosphate.n))

plot(HSD.test(aov(data1$phosphate~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Phosphate")
anova(aov(data1$phosphate~data1$stream+data1$season.yr))

ggplot(data=xphosphate,aes(x=season.yr,y=phosphate.mean))+geom_bar(stat="identity")


ggplot(data=xphosphate, 
       aes(x=season.yr, y=phosphate.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=phosphate.mean, ymax=phosphate.mean+phosphate.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Phosphate ('*'SRP'~~ mg~L^-1*')'))+
  ylim(0,0.04) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=0.0255, label="B", size=3, hjust=0) +
  annotate("text", x=1.95, y=0.023, label="B", size=3, hjust=0)+
  annotate("text", x=2.95, y=0.0355, label="A", size=3, hjust=0)+
  annotate("text", x=2.8, y=0.04, label="P=1.7e-7", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.phosphate1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

################   DIN

xdin.out=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(din.out.mean = mean(din.out, na.rm = TRUE), # na.rm = TRUE to remove missing values
            din.out.sd=sd(din.out, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            din.out.n = sum(!is.na(din.out)), # of observations, excluding NAs. 
            din.out.se=din.out.sd/sqrt(din.out.n))

plot(HSD.test(aov(data1$din.out~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="din.out")
anova(aov(data1$din.out~data1$stream+data1$season.yr))

ggplot(data=xdin.out,aes(x=season.yr,y=din.out.mean))+geom_bar(stat="identity")


ggplot(data=xdin.out, 
       aes(x=season.yr, y=din.out.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=din.out.mean, ymax=din.out.mean+din.out.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Nitrogen ('*'DIN Relative'~~mg~L^-1*')'))+
  ylim(0,0.02) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=0.0115, label="B", size=3, hjust=0) +
  annotate("text", x=1.95, y=0.0085, label="B", size=3, hjust=0)+
  annotate("text", x=2.95, y=0.0175, label="A", size=3, hjust=0)+
  annotate("text", x=2.8, y=0.02, label="P=1.7e-7", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.din1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


################   canopy

xcanopy=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(canopy.mean = mean(canopy, na.rm = TRUE), # na.rm = TRUE to remove missing values
            canopy.sd=sd(canopy, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            canopy.n = sum(!is.na(canopy)), # of observations, excluding NAs. 
            canopy.se=canopy.sd/sqrt(canopy.n))

plot(HSD.test(aov(data1$canopy~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="canopy")
anova(aov(data1$canopy~data1$stream+data1$season.yr))

ggplot(data=xcanopy,aes(x=season.yr,y=canopy.mean))+geom_bar(stat="identity")


ggplot(data=xcanopy, 
       aes(x=season.yr, y=canopy.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=canopy.mean, ymax=canopy.mean+canopy.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab("Canopy (% open)")+
  ylim(0,50) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=34, label="A", size=3, hjust=0) +
  annotate("text", x=1.95, y=44.5, label="A", size=3, hjust=0)+
  annotate("text", x=2.95, y=32., label="A", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.canopy1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


#############  PAR
xpar=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(par.integrative.mean = mean(par.integrative, na.rm = TRUE), # na.rm = TRUE to remove missing values
            par.integrative.sd=sd(par.integrative, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            par.integrative.n = sum(!is.na(par.integrative)), # of observations, excluding NAs. 
            par.integrative.se=par.integrative.sd/sqrt(par.integrative.n))

plot(HSD.test(aov(data1$par.integrative~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="PAR")
anova(aov(data1$par.integrative~data1$stream+data1$season.yr))


ggplot(data=xpar, 
       aes(x=season.yr, y=par.integrative.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=par.integrative.mean, ymax=par.integrative.mean+par.integrative.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('PAR ('*'mol'~~ m^-2~d^-1*')'))+
  ylim(0,2) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=1.8, label="A", size=3, hjust=0) +
  annotate("text", x=1.95, y=.3, label="B", size=3, hjust=0)+
  annotate("text", x=2.95, y=1.35, label="A", size=3, hjust=0)+
  annotate("text", x=2.9, y=2, label="P=0.0010", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.light1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


plot(HSD.test(aov(data1$temp.min~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Min Temp")
anova(aov(data1$temp.min~data1$stream+data1$season.yr))

#############  temp.min
xtemp.min=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(temp.min.mean = mean(temp.min, na.rm = TRUE), # na.rm = TRUE to remove missing values
            temp.min.sd=sd(temp.min, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            temp.min.n = sum(!is.na(temp.min)), # of observations, excluding NAs. 
            temp.min.se=temp.min.sd/sqrt(temp.min.n))

plot(HSD.test(aov(data1$temp.min~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="Min Temp")
anova(aov(data1$temp.min~data1$stream+data1$season.yr))

ggplot(data=xtemp.min, 
       aes(x=season.yr, y=temp.min.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=temp.min.mean, ymax=temp.min.mean+temp.min.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab("Stream Temp. (Minimum °C)")+
  ylim(0,11) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=10.35, label="A", size=3, hjust=0) +
  annotate("text", x=1.95, y=3.3, label="C", size=3, hjust=0)+
  annotate("text", x=2.95, y=8.35, label="B", size=3, hjust=0)+
  annotate("text", x=2.9, y=11, label="P=2.6e-11", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.temp.min1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


#############  cut.mass.m
xcut.mass.m=data1cut%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(cut.mass.m.mean=mean(cut.mass.m,na.rm=TRUE),
            cut.mass.m.max=max(cut.mass.m,na.rm=TRUE),
            cut.mass.m.min=min(cut.mass.m,na.rm=TRUE)) 
            
             # cut.mass.m.sd=sd(cut.mass.m, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            #cut.mass.m.se=cut.mass.m.sd/sqrt(cut.mass.m.n))

ggplot(data=xcut.mass.m,aes(x=season.yr,y=cut.mass.m.mean))+geom_bar(stat="identity")+
geom_errorbar(aes(ymin=1, ymax=13), width=0.2, 
              position=position_dodge(0.9)) 
  

ggplot(data=xcut.mass.m, 
       aes(x=season.yr, y=cut.mass.m.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=cut.mass.m.min, ymax=cut.mass.m.max), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1st Summer '17" = "Summer '17","3rd Summer '18" = "Summer '18"))+
  ylab("Trout Biomass (g/m)")+ ## fix to g-1
  ylim(0,14) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8))

ggsave('N:/Thesis/Rplot.cut.mass.m1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")



#############  gpp.lit
xgpp.lit=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(gpp.lit.mean = mean(gpp.lit, na.rm = TRUE), # na.rm = TRUE to remove missing values
            gpp.lit.sd=sd(gpp.lit, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            gpp.lit.n = sum(!is.na(gpp.lit)), # of observations, excluding NAs. 
            gpp.lit.se=gpp.lit.sd/sqrt(gpp.lit.n))

plot(HSD.test(aov(data1$gpp.lit~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="gpp.lit")
anova(aov(data1$gpp.lit~data1$stream+data1$season.yr))

ggplot(data=xgpp.lit, 
       aes(x=season.yr, y=gpp.lit.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=gpp.lit.mean, ymax=gpp.lit.mean+gpp.lit.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('GPP ('*~O[2]~ m^-2~d^-1*')'))+
  ylim(0,.4) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=.37, label="A", size=3, hjust=0) +
  annotate("text", x=1.95, y=.178, label="B", size=3, hjust=0)+
  annotate("text", x=2.95, y=.2, label="B", size=3, hjust=0)+
  annotate("text", x=2.9, y=.4, label="p= 0.020", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.gpp.lit1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############  er.lit
xer.lit=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(er.lit.mean = abs(mean(er.lit, na.rm = TRUE)), # na.rm = TRUE to remove missing values
            er.lit.sd=sd(er.lit, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            er.lit.n = sum(!is.na(er.lit)), # of observations, excluding NAs. 
            er.lit.se=er.lit.sd/sqrt(er.lit.n))

plot(HSD.test(aov(abs(data1$er.lit)~data1$stream+data1$season.yr), 'data1$season.yr'),ylim=c(0,30),ylab="er.lit")
anova(aov(data1$er.lit~data1$stream+data1$season.yr))

ggplot(data=xer.lit, 
       aes(x=season.yr, y=er.lit.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=er.lit.mean, ymax=er.lit.mean+er.lit.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('ER (Absolute Value'*~O[2]~ m^-2~d^-1*')'))+
  ylim(0,15) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=14.2, label="A", size=3, hjust=0) +
  annotate("text", x=1.95, y=11, label="B", size=3, hjust=0)+
  annotate("text", x=2.91, y=11.9, label="AB", size=3, hjust=0)+
  annotate("text", x=2.9, y=15, label="p= 0.052", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.er.lit1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


#############  pr.lit
xpr.lit=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(pr.lit.mean = mean(pr.lit, na.rm = TRUE), # na.rm = TRUE to remove missing values
            pr.lit.sd=sd(pr.lit, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            pr.lit.n = sum(!is.na(pr.lit)), # of observations, excluding NAs. 
            pr.lit.se=pr.lit.sd/sqrt(pr.lit.n))

plot(HSD.test(aov(data1$pr.lit~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="pr.lit")
anova(aov(data1$pr.lit~data1$stream+data1$season.yr))

ggplot(data=xpr.lit, aes(x=season.yr, y=pr.lit.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=pr.lit.mean, ymax=pr.lit.mean+pr.lit.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab("PR Ratio")+
  ylim(0,.031) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=.031, label="A", size=3, hjust=0) +
  annotate("text", x=1.95, y=.018, label="A", size=3, hjust=0)+
  annotate("text", x=2.95, y=.0215, label="A", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.pr.lit1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############  t.gpp.lit
xt.gpp.lit=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(t.gpp.lit.mean = mean(t.gpp.lit, na.rm = TRUE), # na.rm = TRUE to remove missing values
            t.gpp.lit.sd=sd(t.gpp.lit, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            t.gpp.lit.n = sum(!is.na(t.gpp.lit)), # of observations, excluding NAs. 
            t.gpp.lit.se=t.gpp.lit.sd/sqrt(t.gpp.lit.n))

plot(HSD.test(aov(data1$t.gpp.lit~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="t.gpp.lit")
anova(aov(data1$t.gpp.lit~data1$stream+data1$season.yr))

ggplot(data=xt.gpp.lit, 
       aes(x=season.yr, y=t.gpp.lit.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=t.gpp.lit.mean, ymax=t.gpp.lit.mean+t.gpp.lit.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab("Transformed GPP")+
  ylim(0,.6) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=.95, y=.53, label="A", size=3, hjust=0) +
  annotate("text", x=1.95, y=.39, label="B", size=3, hjust=0)+
  annotate("text", x=2.95, y=.41, label="B", size=3, hjust=0)+
  annotate("text", x=2.4, y=.6, label="Model p< 0.0001", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.t.gpp.lit1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############   t.gpp.depth
ggplot(data=data1, aes(x=depth, y=t.gpp.lit)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab ("Depth (m)")+
  ylab ("Transformed GPP")+
  #ylim(0,0.8)+
  xlim(.04,.13)+
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.title = element_text(size = 6),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=0.04, y=.75,  size=3, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.13"))))+
  annotate("text", x=0.04, y=.71, label="Model p< 0.0001", size=3, hjust=0)


ggsave('N:/Thesis/Rplot.depth.t.gpp1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

##########  t.er.depth
ggplot(data=data1, aes(x=depth, y=t.er.lit)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab ("Depth (m)")+
  ylab ("Transformed ER")+
  #ylim(0,0.8)+
  xlim(.04,.13)+
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.title = element_text(size = 6),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=0.04, y=3.2, size=3, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.36"))))+
  annotate("text", x=0.04, y=3.1, label="Model p< 0.0001", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.depth.t.er1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

##########  t.er.slope
ggplot(data=data1, aes(x=slope, y=t.er.lit)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab ("Slope (%)")+
  ylab ("Transformed ER")+
  #ylim(0,0.8)+
  xlim(2,10)+
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.title = element_text(size = 6),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=2, y=3.2,size=3, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.57"))))+
  annotate("text", x=2, y=3.1, label="Model p< 0.0001", size=3, hjust=0)

ggsave('N:/Thesis/Rplot.slope.t.er1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


##########  t.cut.mass.m.width
ggplot(data=data1, aes(x=width, y=t.cut.mass.m)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab("width (m)")+
  ylab ("Transformed Trout Biomass")+
  #ylim(0,3)+
  xlim(.3,2.2)+
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.title = element_text(size = 6),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=.3, y=2.6, label=("R2= 0.27"), size=3, hjust=0)+
  annotate("text", x=.3, y=2.45, label=("Model P=0.0007"), size=3, hjust=0)
  
  ggsave('N:/Thesis/Rplot.t.cut.mass.m.width1.tiff',
         units="in",
         width=3.25,
         height=3.25,
         dpi=1200,
         compression="lzw")
  
  
##########  t.cut.mass.m.temp.min


ggplot(data=data1, aes(x=temp.min, y=t.cut.mass.m)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab("Stream Temp. (Daily Min. °C)")+
  ylab ("Transformed Trout Biomass")+
  #ylim(0,3)+
  xlim(6.7,10.9)+
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.title = element_text(size = 6),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=9.4, y=2.7, label=("R2=0.19"), size=3, hjust=0)+
  annotate("text", x=9.4, y=2.55, label=("Model P<0.0001"), size=3, hjust=0)

ggsave('N:/Thesis/Rplot.t.cut.mass.m.temp.min1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


##########  t.cut.mass.m.canopy


ggplot(data=data1, aes(x=canopy, y=t.cut.mass.m)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab("Canopy (% open)")+
  ylab ("Transformed Trout Biomass")+
  #ylim(0,3)+
  xlim(4,58)+
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.title = element_text(size = 6),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  #annotate("text", x=9.4, y=2.7, label=("R2=0.19"), size=3, hjust=0)+
  #annotate("text", x=9.4, y=2.55, label=("Model P< 0.0001"), size=3, hjust=0)

ggsave('N:/Thesis/Rplot.t.cut.mass.m.canopy1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")



###########  Jason's fish pop 
ggplot(data=data1, aes(stream, cut.pop.dens, fill=yr))+
geom_bar(stat="identity", position='dodge', color = "black")+  
  geom_errorbar(aes(ymin=cut.pop.dens, ymax=cut.pop.dens+cut.pop.dens.se), width=0.2, position=position_dodge(0.9))+ 
  scale_fill_manual(values=c("gray15","gray65"),name="",
                    breaks=c("A17", "B18"),
                    labels=c("2017", "2018"))+ 
  xlab("Stream")+ 
  scale_x_discrete(limits=c("first","frost","hurley","hovey","blue","swauk","iron","jack","miller","standup"),labels=c("blue" = "Blue","first"="First","frost"="Frost","hovey"="Hovey","hurley"="Hurley","iron"="Iron","jack"="Jack","miller"="Miller","standup"="Standup","swauk"="Swauk"))+
  ylab("Trout Density (fish/m)")+
  ylim(0,2)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
       legend.title=element_text(size=6),
       legend.key=element_blank(),
       legend.position=c(.2,.95),
       legend.text=element_text(size=8),
       legend.background=element_blank(),
       legend.direction="horizontal",
       legend.key.size=unit(0.3, "cm"),
       axis.title.y=element_text(size=8),
       axis.title.x=element_text(size=8),
      axis.text.x=element_text(size=8,angle = 40, hjust = 1))

ggsave('N:/Thesis/Rplot.trout.density1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

###########  Jason's fish mass 
ggplot(data=data1, aes(stream, cut.mass, fill=yr))+
  geom_bar(stat="identity", position='dodge', color = "black")+  
  geom_errorbar(aes(ymin=cut.mass, ymax=cut.mass+cut.mass.se), width=0.2, position=position_dodge(0.9))+ 
  scale_fill_manual(values=c("gray15","gray65"),name="",
                    breaks=c("A17", "B18"),
                    labels=c("2017", "2018"))+ 
  xlab("Stream")+ 
  scale_x_discrete(limits=c("first","frost","hurley","hovey","blue","swauk","iron","jack","miller","standup"),
                   labels=c("blue" = "Blue","first"="First","frost"="Frost","hovey"="Hovey","hurley"="Hurley","iron"="Iron","jack"="Jack","miller"="Miller","standup"="Standup","swauk"="Swauk"))+
  ylab("Mean Trout Mass (g)")+
  ylim(0,41)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(.2,.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8,angle = 40, hjust = 1))

ggsave('N:/Thesis/Rplot.trout.mass1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

###########  Jason's fish mass per meter 
ggplot(data=data1, aes(stream, cut.mass.m, fill=yr))+
  geom_bar(stat="identity", position='dodge', color = "black")+  
  geom_errorbar(aes(ymin=cut.mass.m-cut.mass.m.er, ymax=cut.mass.m+cut.mass.m.er), width=0.2, position=position_dodge(0.9))+ 
  scale_fill_manual(values=c("gray15","gray65"),name="",
                    breaks=c("A17", "B18"),
                    labels=c("2017", "2018"))+ 
  xlab("Stream")+ 
  scale_x_discrete(limits=c("first","frost","hurley","hovey","blue","swauk","iron","jack","miller","standup"),
                   labels=c("blue" = "Blue","first"="First","frost"="Frost","hovey"="Hovey","hurley"="Hurley","iron"="Iron","jack"="Jack","miller"="Miller","standup"="Standup","swauk"="Swauk"))+
  ylab("Trout Biomass (g/m)")+
  ylim(0,18)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(.2,.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8,angle = 40, hjust = 1))

ggsave('N:/Thesis/Rplot.trout.mass.m1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

########## Complete INTERACTION??    

############# t.cut.mass.m.temp.min.canopy 
ggplot(data=dataint, aes(temp.min.cat, t.cut.mass.m.mean, fill=canopy.cat))+
  geom_bar(stat="identity", position='dodge', color = "black")+  
  geom_errorbar(aes(ymin=t.cut.mass.m.min, ymax=t.cut.mass.m.max), width=0.2, position=position_dodge(0.9))+ 
  scale_fill_manual(values=c("gray15","gray65"),name="",
                    breaks=c("c.closed", "c.open"),
                   labels=c("Less open    ", "More Open"))+ 
  xlab("Temperature Category (Minimum °C)")+ 
  scale_x_discrete(limits=c("te.lo","te.mid", "te.hi"),
                   labels=c("te.lo"="Low","te.mid"="Medium", "te.hi"="High"))+
  ylab("Transformed Trout Biomass")+
  ylim(0,2.7)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(.66,.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8,angle = 0, hjust = 1))

ggsave('N:/Thesis/Rplot.interaction1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


###########  Jason's fish mass per meter2 
ggplot(data=data1, aes(stream, cut.mass.m2, fill=yr))+
  geom_bar(stat="identity", position='dodge', color = "black")+  
  geom_errorbar(aes(ymin=cut.mass.m2, ymax=cut.mass.m2+cut.mass.m2.er), width=0.2, position=position_dodge(0.9))+ 
  scale_fill_manual(values=c("gray15","gray65"),name="",
                    breaks=c("A17", "B18"),
                    labels=c("2017", "2018"))+ 
  xlab("Stream")+ 
  scale_x_discrete(limits=c("first","frost","hurley","hovey","blue","swauk","iron","jack","miller","standup"),
                   labels=c("blue" = "Blue","first"="First","frost"="Frost","hovey"="Hovey","hurley"="Hurley","iron"="Iron","jack"="Jack","miller"="Miller","standup"="Standup","swauk"="Swauk"))+
  ylab("Trout Biomass (g/m2)")+
  ylim(0,10)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(.8,.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8,angle = 40, hjust = 1))

ggsave('N:/Thesis/Rplot.biomasssquared1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


###########  Jason's fish pop removed from OnThesis tab
ggplot(data=data1, aes(stream, cut.pop.dens, fill=yr))+
  geom_bar(stat="identity", position='dodge', color = "black")+  
  geom_errorbar(aes(ymin=cut.pop.dens, ymax=cut.pop.dens+cut.pop.dens.se), width=0.2, position=position_dodge(0.9))+ 
  scale_fill_manual(values=c("gray15","gray65"),name="",
                    breaks=c("A17", "B18"),
                    labels=c("2017", "2018"))+ 
  xlab("Stream")+ 
  scale_x_discrete(limits=c("first","frost","hurley","hovey","blue","swauk","iron","jack","miller","standup"),labels=c("blue" = "Blue","first"="First","frost"="Frost","hovey"="Hovey","hurley"="Hurley","iron"="Iron","jack"="Jack","miller"="Miller","standup"="Standup","swauk"="Swauk"))+
  ylab(bquote('Trout Population ('*fish~m^-1*')'))+
  ylim(0,2)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(.2,.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8,angle = 40, hjust = 1))

ggsave('N:/Thesis/Rplot7.trout.density1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")











