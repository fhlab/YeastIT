
library(ggplot2)
library(scales)
library(dplyr)
library(readr)
library(Hmisc)
library(colorscape)
library(viridis)
library(grid)
library(gridExtra)
library(ggpubr)


titv <- read_csv("MutagenicSpectra-titv-withICE_updatedafterreanalysis.csv")
sixtypes <- read_csv("MutagenicSpectra-6types-withICE_updatedafterreanalysis.csv")
rates <- read_csv("MutagenicSpectra-rates-withICE-withPCR_updated.csv")

#transition vs transversion analysis - all methods
p17 <- ggplot(titv, aes(fill=Mutation, y=`%`, x=Method)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#68AD7D", "#A77EAD")) +
  theme_bw(base_size=12) +
  xlab("") + ylab("Mutations (%)") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
   theme(legend.title = element_blank()) +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA)) +
  scale_x_discrete(limits=c("epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "EvolVR", "MutaT7", "eMutaT7", "PACE", "OrthoRep", "TRIDENT", "ICE", "YeastIT1", "YeastIT2"))
p17 <- p17 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p17 <- p17 + theme(legend.position="top") + guides(colour=guide_legend(nrow =1))
p17

#mutation type analysis - all methods
p18 <- ggplot(sixtypes, aes(fill=factor(`Mutation type`, levels=c("AT>GC","GC>AT","AT>CG", "AT>TA", "GC>TA", "GC>CG")), y=`%`, x=Method)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#68AD7D", "#97D2A9", "#A77EAD", "#7c6f91ff", "#8C838D",  "#AE617B")) +
  theme_bw(base_size=12) +
  xlab("") + ylab("Mutations (%)") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA)) 
  p18 <- p18 + scale_x_discrete(limits=c("epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "EvolVR", "MutaT7", "eMutaT7", "PACE", "OrthoRep", "TRIDENT", "ICE", "YeastIT", "YeastIT2"))
  p18 <- p18 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p18 <- p18 + theme(legend.position="top") + guides(colour=guide_legend(nrow =1))
p18


#rates$`Mutation rate` <- as.factor(rates$`Mutation rate`)
p19 <- ggplot(rates, aes(x=Method, y=rates$`Mutation rate`, fill=c("#68AD7D")))
p19 <- p19 + geom_bar(stat="identity") + scale_fill_manual(values = c("#68AD7D")) + theme_bw(base_size=12)
p19 <- p19 + xlab("") + ylab(expression("Mutation rate"~(kb^{-1}~day^{-1}))) 
p19 <- p19 + theme(legend.title = element_blank()) 
p19 <- p19 + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
p19 <- p19 + scale_x_discrete(limits=c("epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "EvolVR", "MutaT7", "eMutaT7", "PACE", "OrthoRep", "TRIDENT", "ICE", "YeastIT1", "YeastIT2"))
p19 <- p19 + theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
p19 <- p19 + theme(legend.position="none") 
p19   
  
p20 <- ggarrange(p17,p18,p19,nrow=3,ncol=1, heights=c(0.8,0.9,1.05))
p20
ggsave("fig2-part-muttypes.svg", width=4.5, height=7.32)


grid.arrange(arrangeGrob(p17, p18, p19, widths=c(1), heights=c(1,1,1), nrow=3))

#old notes below

p5 <- p5 + ylim(-5.5,-0.1) + scale_y_continuous(breaks=c(-6,-5,-4,-3,-2,-1,0)) + facet_wrap(as.factor(growthcounts3$strain)) + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p5 <- p5 + ylab("Log10(Survival fraction)") + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
p5 <- p5 + scale_x_discrete(limits=c("Mut-", "evoAPOBEC1-TadA-T7pol", "pmCDA1-TadA-T7pol", "evoAPOBEC1-T7pol", "pmCDA1-T7pol"))
p5 <- p5 + theme(legend.position = "none") + theme(legend.title = element_blank())
p5 <- p5 +  scale_fill_manual(values = c("#7FC97F", "#a5348f9e", "#7c6f91ff", "#9a9a9ac8",  "#bba137c8",  "#5acf89c8"))
p5
ggsave("mut6poster.png", dpi=600, dev='png', height=4.2, width=6, units="in", bg="transparent")

