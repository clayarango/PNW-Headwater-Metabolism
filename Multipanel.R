#example code for how to make a three panel vertically stacked graph
#each plot is labelled p.1, p.2, p.3
#for each plot, substitute data=x for your data and change axis labels
#reduce to 2 plots for 2 panel figure

#############################################################################################
# Discharge, Canopy, PAR ###########################################################################################
#begin plot 1
pd=position_dodge(0.1)

p.1 = 
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
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
annotate("text", x=.45, y=65, label="A", size=4, hjust=0)

#begin plot 2
pd=position_dodge(0.1)

p.2 = 
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
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
annotate("text", x=.45, y=78, label="B", size=4, hjust=0)

#begin plot 3

pd=position_dodge(0.1)

p.3 = 
  ########################## PAR over time
  ggplot(data=data1, aes(x=season.yr, y=par.integrative)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('PAR ('*mol~photons~m^-2~d^-1*')'))+
  ylim(0,3.9) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9))+
  annotate("text", x=0.955, y=3.73, label="A", size=3.5, hjust=0) +
  annotate("text", x=1.955, y=.78, label="B", size=3.5, hjust=0)+
  annotate("text", x=2.955, y=2.39, label="A", size=3.5, hjust=0)+
  annotate("text", x=2.950, y=3.73, label="p=0.001", size=3.5, hjust=0)+
  annotate("text", x=.45, y=3.73, label="C", size=4, hjust=0)
# P-value from tukey hsd

#Make First three panel figure
gA <- ggplotGrob(p.1)  # set up figure
gB <- ggplotGrob(p.2)  # set up figure
gC <- ggplotGrob(p.3)

#ensure widths line up properly
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5], gC$widths[2:5])  # set up figure
gA$widths[2:5] <- as.list(maxWidth)  # set up figure
gB$widths[2:5] <- as.list(maxWidth)  # set up figure
gC$widths[2:5] <- as.list(maxWidth)

#save plot (need to update filename = )
tiff(filename ='N:/Thesis/Rplot1.discharge.canopy.par.tiff', #open plotting device
     width = 4.0,
     height = 6.5,
     units = "in",
     res = 1200,
     compression = "lzw")
grid.arrange(gA, gB, gC, nrow=3, ncol=1)  # push plot to device
dev.off()  # close device

#################################################################################################
# DIN, SRP, DOC ##################################################################################

#begin plot 1
pd=position_dodge(0.1)

p.1 = 
  ########################## DIN over time
  ggplot(data=data1, aes(x=season.yr, y=din.out)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Relative Nitrogen ('*DIN~mg~N~L^-1*')'))+
  #ylim(0,.0185)+ 
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  annotate("text", x=0.955, y=.0164, label="B", size=3.5, hjust=0) +
  annotate("text", x=1.955, y=.0158, label="B", size=3.5, hjust=0)+
  annotate("text", x=2.955, y=.0183, label="A", size=3.5, hjust=0)+
  annotate("text", x=1.2, y=.0183, label="p=4.6e-7", size=3.5, hjust=0)+
  annotate("text", x=.45, y=0.0183, label="A", size=4, hjust=0)
# p value from tukey hsd
#begin plot 2
pd=position_dodge(0.1)

p.2 = 
  ########################## SRP over time
  ggplot(data=data1, aes(x=season.yr, y=phosphate)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Phosphate ('*SRP~mg~P~L^-1*')'))+
  #ylim(0,.064) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  annotate("text", x=0.955, y=.0455, label="B", size=3.5, hjust=0) +
  annotate("text", x=1.955, y=.0475, label="B", size=3.5, hjust=0)+
  annotate("text", x=2.955, y=.0631, label="A", size=3.5, hjust=0)+
  annotate("text", x=1.2, y=.0631, label="p=1.7e-7", size=3.5, hjust=0)+
annotate("text", x=.45, y=0.062, label="B", size=4, hjust=0)
# p value from tukey hsd

#begin plot 3

pd=position_dodge(0.1)

p.3 = 
  ########################## DOC over time
  ggplot(data=data1, aes(x=season.yr, y=carbon))+  
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('Carbon ('*DOC~mg~C~L^-1*')'))+
  #ylim(0,14.3) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9))+
  annotate("text", x=0.955, y=7.5, label="B", size=3.5, hjust=0) +
  annotate("text", x=1.955, y=9.6, label="B", size=3.5, hjust=0)+
  annotate("text", x=2.955, y=13.8, label="A", size=3.5, hjust=0)+
  annotate("text", x=1.2, y=13.8, label="p=1.2e-4", size=3.5, hjust=0)+
  annotate("text", x=.45, y=13.6, label="C", size=4, hjust=0)
# p value from tukey hsd
  
#Make First three panel figure
gA <- ggplotGrob(p.1)  # set up figure
gB <- ggplotGrob(p.2)  # set up figure
gC <- ggplotGrob(p.3)

#ensure widths line up properly
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5], gC$widths[2:5])  # set up figure
gA$widths[2:5] <- as.list(maxWidth)  # set up figure
gB$widths[2:5] <- as.list(maxWidth)  # set up figure
gC$widths[2:5] <- as.list(maxWidth)

#save plot (need to update filename = )
tiff(filename ='N:/Thesis/Rplot2.din.srp.doc.tiff', #open plotting device
     width = 4.0,
     height = 6.5,
     units = "in",
     res = 1200,
     compression = "lzw")
grid.arrange(gA, gB, gC, nrow=3, ncol=1)  # push plot to device
dev.off()  # close device

#############################################################################################
# GPP #########################################################################################

ggplot(data=data1, aes(x=season.yr, y=gpp.lit)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('GPP ('*g~O[2]~  m^-2~d^-1*')'))+
  #ylim(0,.4) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9))+
  annotate("text", x=.94, y=.725, label="A", size=3.5, hjust=0) +
  annotate("text", x=1.95, y=.30, label="B", size=3.5, hjust=0)+
  annotate("text", x=2.95, y=.37, label="B", size=3.5, hjust=0)+
  annotate("text", x=2.7, y=.725, label="p<0.0001", size=3.5, hjust=0)
# P value from LME

ggsave('N:/Thesis/Rplot3.gpp.lit.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############################################################################################
# t.GPP.depth #########################################################################################
ggplot(data=data1, aes(x=depth, y=t.gpp.lit)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab ("Depth (m)")+
  ylab (expression(sqrt(ln(GPP+1))))+
  #ylim(0,0.8)+
  xlim(.04,.13)+
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9),
        legend.title = element_text(size = 9),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=9),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=0.04, y=.73,  size=3.5, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.13"))))+
  annotate("text", x=0.04, y=.68, label="p<0.0001", size=3.5, hjust=0)
# p value from LME

ggsave('N:/Thesis/Rplot4.depth.t.gpp.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############################################################################################
# ER #########################################################################################

ggplot(data=data1, aes(x=season.yr, y=abs(er.lit))) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65")  +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Sampling Period") +
  scale_x_discrete(labels=c("1sum17" = "Summer '17", "2fall17" = "Fall '17","3sum18" = "Summer '18"))+
  ylab(bquote('|ER| ('*g~O[2]~ m^-2~d^-1*')'))+
  #ylim(0,15) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9))
  #annotate("text", x=.95, y=24.9, label="A", size=3.5, hjust=0) +
  #annotate("text", x=1.95, y=14.9, label="B", size=3.5, hjust=0)+
  #annotate("text", x=2.91, y=21.85, label="AB", size=3.5, hjust=0)+
  #annotate("text", x=2.7, y=24.9, label="p=0.052", size=3.5, hjust=0)
# P value from tukey hsd

ggsave('N:/Thesis/Rplot5.er.lit.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

################################################################################################
##########  t.er.depth ######################################################################
ggplot(data=data1, aes(x=depth, y=t.er.lit)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab ("Depth (m)")+
  ylab (bquote('ln (|'*ER*'|+1)'))+
  #ylim(0,0.8)+
  xlim(.04,.13)+
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9),
        legend.title = element_text(size = 9),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=9),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=0.04, y=3.2, size=3.5, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.36"))))+
  annotate("text", x=0.04, y=3.1, label="p<0.0001", size=3.5, hjust=0)
# p value from LME

ggsave('N:/Thesis/Rplot6.t.er.depth.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")
###############################################################################################
##########  t.er.slope ##########################################################################
ggplot(data=data1, aes(x=slope, y=t.er.lit)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab ("Slope (%)")+
  ylab (bquote('ln (|'*ER*'|+1)'))+
  #ylim(0,0.8)+
  xlim(2,10)+
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9),
        legend.title = element_text(size = 9),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=9),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=2, y=3.2,size=3.5, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.57"))))+
  annotate("text", x=2, y=3.1, label="p<0.0001", size=3.5, hjust=0)
#p value from lme

ggsave('N:/Thesis/Rplot7.t.er.slope.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############################################################################################
# gpp.er #########################################################################################

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
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9),
        legend.title = element_text(size = 9),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=9),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=0, y=24,size=3.5, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.41"))))+
  annotate("text", x=0, y=22.5, label="p= 0.00026", size=3.5, hjust=0)
# p value from ?

ggsave('N:/Thesis/Rplot8.gpp.er.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############################################################################################
# fish #########################################################################################

#begin plot 1
pd=position_dodge(0.1)

p.1 = 
  ###########  Trout pop density
  ggplot(data=data1, aes(stream, cut.pop.dens, fill=yr))+
  geom_bar(stat="identity", position='dodge', color = "black")+  
  geom_errorbar(aes(ymin=cut.pop.dens, ymax=cut.pop.dens+cut.pop.dens.se), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray15","gray65"),name="",
                    breaks=c("A17", "B18"),
                    labels=c("2017", "2018"))+ 
  xlab("Stream")+ 
  scale_x_discrete(limits=c("first","frost","hurley","hovey","blue","swauk","iron","jack","miller","standup"),
                   labels=c("blue" = "Blue","first"="First","frost"="Frost","hovey"="Hovey","hurley"="Hurley","iron"="Iron","jack"="Jack","miller"="Miller","standup"="Standup","swauk"="Swauk"))+
  ylab(bquote('Trout Population ('*fish~m^-1*')'))+
  #ylim(0,41)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(.27,.94),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
annotate("text", x=.50, y=1.84, label="A", size=4, hjust=0)

#begin plot 2
pd=position_dodge(0.1)

p.2 = 
  ########### indiv trout mass 
  ggplot(data=data1, aes(stream, cut.mass, fill=yr))+
  geom_bar(stat="identity", position='dodge', color = "black")+  
  geom_errorbar(aes(ymin=cut.mass, ymax=cut.mass+cut.mass.se), width=0.2, position=position_dodge(0.9))+ 
  scale_fill_manual(values=c("gray15","gray65"),name="",
                    breaks=c("A17", "B18"),
                    labels=c("2017", "2018"))+ 
  xlab("Stream")+ 
  scale_x_discrete(limits=c("first","frost","hurley","hovey","blue","swauk","iron","jack","miller","standup"),
                   labels=c("blue" = "Blue","first"="First","frost"="Frost","hovey"="Hovey","hurley"="Hurley","iron"="Iron","jack"="Jack","miller"="Miller","standup"="Standup","swauk"="Swauk"))+
  ylab("Individual Trout Mass (g)")+
  #ylim(0,41)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
annotate("text", x=.50, y=39.5, label="B", size=4, hjust=0)

#begin plot 3

pd=position_dodge(0.1)

p.3 = 
  ###########  fish mass per meter2 
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
  #ylim(0,10)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9,angle =30, hjust = 1))+
  annotate("text", x=.50, y=9.4, label="C", size=4, hjust=0)

#Make First three panel figure
gA <- ggplotGrob(p.1)  # set up figure
gB <- ggplotGrob(p.2)  # set up figure
gC <- ggplotGrob(p.3)

#ensure widths line up properly
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5], gC$widths[2:5])  # set up figure
gA$widths[2:5] <- as.list(maxWidth)  # set up figure
gB$widths[2:5] <- as.list(maxWidth)  # set up figure
gC$widths[2:5] <- as.list(maxWidth)

#save plot (need to update filename = )
tiff(filename ='N:/Thesis/Rplot9.fish.tiff', #open plotting device
     width = 4.0,
     height = 6.5,
     units = "in",
     res = 1200,
     compression = "lzw")
grid.arrange(gA, gB, gC, nrow=3, ncol=1)  # push plot to device

dev.off()  # close device

#############################################################################################
# fish.catchment #########################################################################################
ggplot(data=data1, aes(x=basin, y=t.cut.mass.m2)) + 
  geom_boxplot(position=position_dodge(),color = "black",fill="gray65") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=2.5)+ 
  xlab("Catchment") +
  scale_x_discrete(limits=c("taneum","swauk","teanaway"),
                   labels=c("swauk" = "Swauk", "taneum" = "Taneum","teanaway" = "Teanaway"))+
  ylab("ln(Trout Biomass+1)")+
  #ylim(0,2.3) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9))+
  annotate("text", x=0.65, y=2.2, label="p=0.0007", size=3.5, hjust=0)+
  annotate("text", x=1.95, y=2.30, label="A", size=3, hjust=0)+
  annotate("text", x=2.90, y=2.10, label="AB", size=3, hjust=0)+
  annotate("text", x=0.95, y=1.86, label="B", size=3, hjust=0)
# p value from GLS

ggsave('N:/Thesis/Rplot10.t.cut.mass.catchment.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############################################################################################
# fish.interaction #########################################################################################
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
  ylab("ln(Trout Biomass+1)")+
  #ylim(0,2.5)+
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=9),
        legend.key=element_blank(),
        legend.position=c(.8,.94),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        legend.direction="vertical",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9,angle = 0, hjust = 1))+
  annotate("text", x=2.6, y=.1, label="p=0.0071", size=3.5, hjust=0)
# p value from GLS interaction

ggsave('N:/Thesis/Rplot11.fish.interaction.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#############################################################################################
# gpp.fish #########################################################################################
ggplot(data=data1, aes(x=gpp.lit, y=cut.mass.m2)) +
  geom_point(size=2, shape=16,color="gray30") +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab(bquote('GPP ('*g~O[2]~ m^-2~d^-1*')'))+
  ylab(bquote('Trout Biomass ('*g~m^-2*')'))+
  ylim(0,8.5)+
  xlim(0,.71)+
  #stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_text(size=9),
        legend.title = element_text(size = 9),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=9),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=0.5, y=8.25,  size=3.5, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.012"))))+
  annotate("text", x=0.5, y=7.65, label="p= 0.28", size=3.5, hjust=0)

ggsave('N:/Thesis/Rplot12.gpp.fish.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")
#############################################################################################
##########################################################################################
setwd("N:/Thesis")
swaukexample=read.csv("N:/Thesis/swauk18example2.csv")




ggplot() +
  geom_point(data=swaukexample, aes(x=seq, y=value,group="oxy", 
                                    color=(oxy),size=(oxy)))+
  scale_color_manual(labels=c(" Modeled  "," Observed   "),values=c("black", "gray50"))+
  
  
  scale_size_manual(values=c(.5,.9))+
  xlab("Time (Midnight to Midnight)")+
  ylab(bquote(''*O[2]~(mg~L^-1*'')))+
#  #ylim(8.99,9.81)+
#  #xlim(0,300)+
  #stat_smooth(method="lm", se=F, col="black") 
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),  
        legend.position = c(0.8, 1),  
        legend.text=element_text(size=9),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm"))+
  annotate("text", x=1, y=9.79,  size=3.5, hjust=0,parse = TRUE, label = as.character(expression(paste(R[adj]^"2","=0.99"))))
#  annotate("text", x=0.5, y=7.65, label="p= 0.28", size=3.5, hjust=0)

ggsave('N:/Thesis/Rplot.oxy.mod.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")







