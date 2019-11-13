install.packages("tidyverse")
install.packages("nlme4")
install.packages("nortest")
install.packages("multcomp")
install.packages("MASS")
install.packages("gridExtra")
install.packages("vegan")
install.packages("agricolae")

library(tidyverse)
library(nlme4)
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
data1$t.cut.mass.m2=log(1+data1$cut.mass.m2)

#############  gpp.lit
plot(HSD.test(aov(data1$gpp.lit~data1$stream+data1$season.yr), 'data1$season.yr'),ylab="gpp.lit")
anova(aov(data1$gpp.lit~data1$stream+data1$season.yr))

xgpp.lit=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(gpp.lit.mean = mean(gpp.lit, na.rm = TRUE), # na.rm = TRUE to remove missing values
            gpp.lit.sd=sd(gpp.lit, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            gpp.lit.n = sum(!is.na(gpp.lit)), # of observations, excluding NAs. 
            gpp.lit.se=gpp.lit.sd/sqrt(gpp.lit.n))

ggplot(data=xgpp.lit, 
       aes(x=season.yr, y=gpp.lit.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=gpp.lit.mean, ymax=gpp.lit.mean+gpp.lit.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('GPP ('*g~O[2]~  m^-2~d^-1*')'))+
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
  annotate("text", x=2.4, y=.4, label="LME p< 0.0001", size=3, hjust=0)

ggsave('N:/Thesis/Rplot1.gpp.lit1.tiff',
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
  ylab (expression(sqrt(ln(1+GPP))))+
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
  annotate("text", x=0.04, y=.71, label="LME p< 0.0001", size=3, hjust=0)

ggsave('N:/Thesis/Rplot2.depth.t.gpp1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############  er.lit
plot(HSD.test(aov(abs(data1$er.lit)~data1$stream+data1$season.yr), 'data1$season.yr'),ylim=c(0,30),ylab="er.lit")
anova(aov(data1$er.lit~data1$stream+data1$season.yr))

xer.lit=data1%>% group_by(season.yr)%>%   # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(er.lit.mean = abs(mean(er.lit, na.rm = TRUE)), # na.rm = TRUE to remove missing values
            er.lit.sd=sd(er.lit, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            er.lit.n = sum(!is.na(er.lit)), # of observations, excluding NAs. 
            er.lit.se=er.lit.sd/sqrt(er.lit.n))

ggplot(data=xer.lit, 
       aes(x=season.yr, y=er.lit.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=er.lit.mean, ymax=er.lit.mean+er.lit.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('|ER| ('*g~O[2]~ m^-2~d^-1*')'))+
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
  annotate("text", x=2.2, y=15, label="Tukey HSD p= 0.052", size=3, hjust=0)

ggsave('N:/Thesis/Rplot3.er.lit1.tiff',
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
  ylab (bquote('ln (1+|'*ER*'|)'))+
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
  annotate("text", x=0.04, y=3.1, label="LME p< 0.0001", size=3, hjust=0)

ggsave('N:/Thesis/Rplot4.depth.t.er1.tiff',
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
  ylab (bquote('ln (1+|'*ER*'|)'))+
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
  annotate("text", x=2, y=3.1, label="LME p< 0.0001", size=3, hjust=0)

ggsave('N:/Thesis/Rplot5.slope.t.er1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

##########  gpp.er
ggplot(data=data1, aes(x=gpp.lit, y=abs(er.lit))) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab(bquote('GPP ('*g~O[2]~ m^-2~d^-1*')'))+
  ylab(bquote('|ER| ('*g~O[2]~ m^-2~d^-1*')'))+
  #ylim(0,0.8)+
  xlim(0,.71)+
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
  annotate("text", x=0, y=24,size=3, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.41"))))+
  annotate("text", x=0, y=22.5, label="p= 0.00026", size=3, hjust=0)

ggsave('N:/Thesis/Rplot6.gpp.er1.tiff',
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
  ylab("Mass of Individual Fish (g)")+
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

ggsave('N:/Thesis/Rplot7.trout.mass1.tiff',
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
  ylab(bquote('Trout Biomass ('*g~m^-2*')'))+
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

ggsave('N:/Thesis/Rplot8.biomasssquared1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############  t.cut.mass.m.basin
ggplot(data=data1, aes(x=basin, y=t.cut.mass.m2)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Catchment") +
  scale_x_discrete(limits=c("taneum","swauk","teanaway"),
                   labels=c("swauk" = "Swauk", "taneum" = "Taneum","teanaway" = "Teanaway"))+
  ylab("ln(Trout Biomass+1)")+
  ylim(0,2.3) +
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
        axis.text.x=element_text(size=8))+
  annotate("text", x=0.6, y=2.3, label="p=0.0007", size=3, hjust=0)+
annotate("text", x=1.95, y=2.30, label="A", size=3, hjust=0)+
annotate("text", x=2.90, y=2.10, label="AB", size=3, hjust=0)+
annotate("text", x=0.95, y=1.86, label="B", size=3, hjust=0)

ggsave('N:/Thesis/Rplot9.t.cut.mass.m.basin1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

############# t.cut.mass.m.temp.min.canopy 
ggplot(data=dataint, aes(temp.min.cat, t.cut.mass.m2, fill=canopy.cat))+
  geom_boxplot(width=0.7, position=position_dodge(0.9)) +
  scale_fill_manual(values=c("gray35","gray75"),name="",
                    breaks=c("c.closed", "c.open"),
                   labels=c("Less open   ", "More open")) +
  xlab("Stream Temperature Category (Daily Min. °C)") +
  scale_x_discrete(limits=c("te.lo","te.mid", "te.hi"),
                   labels=c("te.lo"="Low","te.mid"="Mid", "te.hi"="Hi"))+
  # step =1.4 degrees between temp categories and 25.3 % openness between canopy categories
  # less open=4.9-30.2 %, more open=30.2-55.5 %, lo=6.8-8.1, mid=8.1-9.5, hi=9.5-10.9 deg
  ylab("ln(1+Trout Biomass)")+
  ylim(0,2.5)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(.83,.9),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="vertical",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8,angle = 0, hjust = 1))+
  annotate("text", x=1.8, y=2.5, label="GLS Interaction p= 0.0071", size=3, hjust=0)

ggsave('N:/Thesis/Rplot10.interaction1.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############   gpp.trout
ggplot(data=data1, aes(x=gpp.lit, y=cut.mass.m2)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab(bquote('GPP ('*g~O[2]~ m^-2~d^-1*')'))+
  ylab(bquote('Trout Biomass ('*g~m^-2*')'))+
  ylim(0,8.5)+
  xlim(0,.71)+
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
  annotate("text", x=0.55, y=8.25,  size=3, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.012"))))+
  annotate("text", x=0.55, y=7.65, label="p= 0.28", size=3, hjust=0)

ggsave('N:/Thesis/Rplot11.gpp.trout.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

########################## Discharge over time
ggplot(data=data1, aes(x=season.yr, y=discharge)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
   xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Stream Discharge ('*L~ s^-1*')'))+
  ylim(0,66) +
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
#  annotate("text", x=2.4, y=2.3, label="GLS p= 0.0007", size=3, hjust=0)

ggsave('N:/Thesis/Rplot12.discharge.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


########################## canopy over time
ggplot(data=data1, aes(x=season.yr, y=canopy)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab("Canopy Openness (%)")+
  ylim(0,79) +
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
#  annotate("text", x=2.4, y=2.3, label="GLS p= 0.0007", size=3, hjust=0)

ggsave('N:/Thesis/Rplot13.canopy.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


########################## PAR over time
ggplot(data=data1, aes(x=season.yr, y=par.integrative)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('PAR ('*mol~photons~m^-2~d^-1*')'))+
  ylim(0,3.7) +
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
        axis.text.x=element_text(size=8))+
annotate("text", x=.956, y=3.64, label="A", size=3, hjust=0) +
  annotate("text", x=1.957, y=.67, label="B", size=3, hjust=0)+
  annotate("text", x=2.958, y=2.36, label="A", size=3, hjust=0)+
  annotate("text", x=2.1, y=3.64, label="Tukey HSD p= 0.001", size=3, hjust=0)

ggsave('N:/Thesis/Rplot14.PAR.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")


########################## DIN over time
ggplot(data=data1, aes(x=season.yr, y=din.out)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Relative Nitrogen ('*DIN~mg~L^-1*')'))+
  ylim(.0020,.0185)+ 
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
        axis.text.x=element_text(size=8))+
  annotate("text", x=.956, y=.0164, label="B", size=3, hjust=0) +
  annotate("text", x=1.957, y=.0158, label="B", size=3, hjust=0)+
  annotate("text", x=2.954, y=.0182, label="A", size=3, hjust=0)+
  annotate("text", x=.65, y=.0182, label="Tukey HSD p= 4.6e-7", size=3, hjust=0)

ggsave('N:/Thesis/Rplot15.DIN.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

########################## SRP over time
ggplot(data=data1, aes(x=season.yr, y=phosphate)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Phosphate ('*SRP~mg~L^-1*')'))+
  ylim(.004,.064) +
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
        axis.text.x=element_text(size=8))+
  annotate("text", x=.956, y=.0455, label="B", size=3, hjust=0) +
  annotate("text", x=1.957, y=.0472, label="B", size=3, hjust=0)+
  annotate("text", x=2.954, y=.0631, label="A", size=3, hjust=0)+
  annotate("text", x=.65, y=.0631, label="Tukey HSD p= 1.7e-7", size=3, hjust=0)

ggsave('N:/Thesis/Rplot16.SRP.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

########################## DOC over time
ggplot(data=data1, aes(x=season.yr, y=carbon))+  
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Carbon ('*DOC~mg~L^-1*')'))+
  ylim(0,14.3) +
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
        axis.text.x=element_text(size=8))+
  annotate("text", x=.956, y=7.3, label="B", size=3, hjust=0) +
  annotate("text", x=1.957, y=9.6, label="B", size=3, hjust=0)+
  annotate("text", x=2.954, y=13.8, label="A", size=3, hjust=0)+
  annotate("text", x=.65, y=13.8, label="Tukey HSD p= 1.2e-4", size=3, hjust=0)

ggsave('N:/Thesis/Rplot17.DOC.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

