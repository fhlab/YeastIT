
library(ggplot2)
library(scales)
library(dplyr)
library(readr)
library(Hmisc)
library(RColorBrewer)
library(grid)
library(gridExtra)

aacategories_stats <- read.csv("stats.csv", header=TRUE)

aacategories_stats$Filename[aacategories_stats$Filename=="MutazymeII"] <- "epPCR (Mutazyme II)"
aacategories_stats$Filename[aacategories_stats$Filename=="Ideal"] <- "Equal Probabilities"
aacategories_stats$Filename[aacategories_stats$Filename=="MutaT7"] <- "(e)MutaT7"
aacategories_stats$Filename[aacategories_stats$Filename=="Taq"] <- "epPCR (Taq, Mn2+)"

aacategories_stats$ANYtoPCoverPCtoNP <- (aacategories_stats$NP.to.P + aacategories_stats$NP.to.B + aacategories_stats$NP.to.A +  aacategories_stats$P.to.A + aacategories_stats$P.to.B + aacategories_stats$A.to.B + aacategories_stats$A.to.P + aacategories_stats$P.to.P +  aacategories_stats$A.to.A + aacategories_stats$B.to.B +  aacategories_stats$B.to.A)    / (aacategories_stats$P.to.NP + aacategories_stats$A.to.NP + aacategories_stats$B.to.NP)
aacategories_stats$ANYtoPCwithoutselfoverPCtoNP <- (aacategories_stats$NP.to.P + aacategories_stats$NP.to.B + aacategories_stats$NP.to.A +  aacategories_stats$P.to.A + aacategories_stats$P.to.B + aacategories_stats$A.to.B + aacategories_stats$A.to.P    +  aacategories_stats$B.to.A)    / (aacategories_stats$P.to.NP + aacategories_stats$A.to.NP + aacategories_stats$B.to.NP)
aacategories_stats$acidoverbasewithoutneutralizing <- (aacategories_stats$B.to.A + aacategories_stats$NP.to.A  +  aacategories_stats$P.to.A + aacategories_stats$A.to.A)    / (aacategories_stats$P.to.B + aacategories_stats$A.to.B + aacategories_stats$NP.to.B + aacategories_stats$B.to.B)
aacategories_stats$acidoverbase <- (aacategories_stats$B.to.A + aacategories_stats$NP.to.A  +  aacategories_stats$P.to.A + aacategories_stats$A.to.A + aacategories_stats$B.to.P + aacategories_stats$B.to.NP)    / (aacategories_stats$P.to.B + aacategories_stats$A.to.B + aacategories_stats$NP.to.B + aacategories_stats$B.to.B + aacategories_stats$A.to.NP + aacategories_stats$A.to.P)


p1 <- ggplot(aacategories_stats, aes(x=reorder(Filename, -ANYtoPCwithoutselfoverPCtoNP), y=ANYtoPCwithoutselfoverPCtoNP, fill=c("#7FB9DA")))
p1 <- p1 + geom_bar(stat="identity") + scale_fill_manual(values = c("#7FB9DA")) + theme_grey(base_size=13)
p1 <- p1 + xlab("") + ylab("Propensity for functional innovation") 
p1 <- p1 + theme(legend.title = element_blank()) 
p1 <- p1 + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
#p1 <- p1 + scale_x_discrete(limits=c("Equal Probabilities", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "MutaT7", "PACE", "OrthoRep", "TRIDENT", "YeastIT1", "YeastIT2"))
#uncomment above for original order
p1 <- p1 + theme(legend.position="none") 
p1 <- p1 + theme(axis.text.x = element_text(angle = 50, hjust = 1))
p1 

p2 <- ggplot(aacategories_stats, aes(x=reorder(Filename, -acidoverbase), y=acidoverbase, fill=c("#7FB9DA")))
p2 <- p2 + geom_bar(stat="identity") + scale_fill_manual(values = c("#7FB9DA")) + theme_grey(base_size=13)
p2 <- p2 + xlab("") + ylab("Propensity for acidity") 
p2 <- p2 + theme(legend.title = element_blank()) 
p2 <- p2 + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
#p2 <- p2 + scale_x_discrete(limits=c("Equal Probabilities", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "MutaT7", "PACE", "OrthoRep", "TRIDENT", "YeastIT1", "YeastIT2"))
#uncomment above for original order
p2 <- p2 + theme(legend.position="none") 
p2 <- p2 + theme(axis.text.x = element_text(angle = 50, hjust = 1))
p2 

grid.arrange(arrangeGrob(p1, p2, widths=c(1,1), nrow=1))

