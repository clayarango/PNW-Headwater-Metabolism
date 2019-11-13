#example code for how to make a three panel vertically stacked graph
#each plot is labelled p.1, p.2, p.3
#for each plot, substitute data=x for your data and change axis labels
#reduce to 2 plots for 2 panel figure

#begin plot 1
pd=position_dodge(0.1)

p.1 = 
  ggplot(data=x, aes(x=cat.time, y=emmean.raw, fill=budworm)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=emmean.raw, ymax=emmean.raw+SE.raw), width=0.2, 
                position=position_dodge(0.9)) + 
  xlab("Sample Event") +
  ylab(expression(Throughfall~VWM~DIN~(mg~N~L^{-1}))) +
  ylim(0,1.0) +
  scale_fill_manual(name="Budworm Activity", values=c("black", "white")) +
  expand_limits(y=0.02) +
  theme_bw() +
  theme(legend.justification=c(0.03,0.6),
        legend.position=c(0.075,0.85),
        legend.title=element_text(size= 12),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  annotate("text", x=1, y=0.25, label="+", size=4) +
  #annotate("text", x=2.1, y=0.25, label="bc", size=4) +
  #annotate("text", x=3.1, y=0.23, label="c", size=4) +
  #annotate("text", x=4, y=0.15, label="a", size=4) +
  annotate("text", x=5, y=0.2, label="+", size=4) +
  #annotate("text", x=6.1, y=0.17, label="b", size=4) +
  #annotate("text", x=6.8, y=0.22, label="c", size=4) +
  #annotate("text", x=8.3, y=0.9, label="c", size=4) +
  annotate("text", x=9, y=0.2, label="+", size=4) +
  annotate("text", x=10, y=0.2, label="+", size=4) +
  annotate("text", x=10, y=0.95, label="A", size=6)
#annotate("text", x=1.25, y=0.85, label="budworm,", size=4, hjust=0) +
#annotate("text", x=1.45, y=0.79, label="p=0.68", size=4, hjust=0) +
#annotate("text", x=1.25, y=0.73, label="sample event,", size=4, hjust=0) +
#annotate("text", x=1.45, y=0.67, label="p=0.003", size=4, hjust=0) +
#annotate("text", x=1.25, y=0.61, label="budworm*sample event,", size=4, hjust=0) +
#annotate("text", x=1.45, y=0.55, label="p=0.010", size=4, hjust=0)



#begin plot 2
pd=position_dodge(0.1)

p.2 = 
  ggplot(data=x, aes(x=cat.time, y=emmean.raw, fill=budworm)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=emmean.raw, ymax=emmean.raw+SE.raw), width=0.2, 
                position=position_dodge(0.9)) +
  xlab("Sample Event") +
  ylab(expression(Throughfall~VWM~PO[4]^{"3-"}~(mg~P~L^{-1}))) +
  scale_fill_manual(name="Budworm Activity", values=c("black", "white")) +
  theme_bw() +
  theme(legend.position="none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  annotate("text", x=1, y=10, label="+", size=4) +
  #annotate("text", x=2.2, y=10, label="bc", size=4) +
  #annotate("text", x=3.2, y=12, label="cd", size=4) +
  #annotate("text", x=4.2, y=14, label="cd", size=4) +
  #annotate("text", x=5.2, y=8, label="a", size=4) +
  #annotate("text", x=6.2, y=8, label="a", size=4) +
  #annotate("text", x=6.8, y=12, label="de", size=4) +
  annotate("text", x=8, y=100, label="+", size=4) +
  annotate("text", x=9, y=8, label="*", size=4) +
  #annotate("text", x=10.1, y=14, label="cd", size=4) +
  annotate("text", x=10, y=95, label="B", size=6)
#annotate("text", x=1.25, y=90, label="budworm,", size=4, hjust=0) +
#annotate("text", x=1.45, y=81, label="p=0.0014", size=4, hjust=0) +
#annotate("text", x=1.25, y=72, label="sample event,", size=4, hjust=0) +
#annotate("text", x=1.45, y=63, label="p<0.0001", size=4, hjust=0) +
#annotate("text", x=1.25, y=54, label="budworm*sample event,", size=4, hjust=0) +
#annotate("text", x=1.45, y=45, label="p=0.0001", size=4, hjust=0) 


#begin plot 3

pd=position_dodge(0.1)

p.3 = 
  ggplot(data=x, aes(x=cat.time, y=emmean, fill=budworm)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=emmean, ymax=emmean+SE), width=0.2, 
                position=position_dodge(0.9)) +
  xlab("Sample Event") +
  #ylim(0,700) +
  ylab(expression(Throughfall~VWM~DOC~(mg~C~L^{-1}))) +
  scale_fill_manual(name="Budworm Activity", values=c("black", "white")) +
  expand_limits(y=0.02) +
  theme_bw() +
  theme(legend.position="none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title=element_text(size= 12),
        legend.text=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  #annotate("text", x=1.2, y=270, label="ab", size=4) +
  annotate("text", x=2, y=300, label="+", size=4) +
  #annotate("text", x=3.2, y=250, label="ab", size=4) +
  annotate("text", x=4, y=200, label="+", size=4) +
  #annotate("text", x=5.2, y=120, label="a", size=4) +
  #annotate("text", x=5.8, y=130, label="ab", size=4) +
  #annotate("text", x=6.8, y=350, label="b", size=4) +
  annotate("text", x=8, y=600, label="*", size=4) +
  #annotate("text", x=9.2, y=200, label="ab", size=4) +
  #annotate("text", x=9.8, y=230, label="a", size=4) +
  annotate("text", x=10, y=570, label="C", size=6)
#annotate("text", x=1.25, y=670, label="budworm,", size=4, hjust=0) +
#annotate("text", x=1.45, y=610, label="p=0.43", size=4, hjust=0) +
#annotate("text", x=1.25, y=550, label="sample event,", size=4, hjust=0) +
#annotate("text", x=1.45, y=490, label="p<0.0001", size=4, hjust=0) +
#annotate("text", x=1.25, y=430, label="budworm*sample event,", size=4, hjust=0) +
#annotate("text", x=1.45, y=370, label="p=0.0002", size=4, hjust=0)

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
tiff(filename = 'figures/new figures/TF_DIN_SRP_DOC_Concentration.tiff', #open plotting device
     width = 6.5,
     height = 9.0,
     units = "in",
     res = 1200,
     compression = "lzw")
grid.arrange(gA, gB, gC, nrow=3, ncol=1)  # push plot to device
dev.off()  # close device