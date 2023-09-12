
library(ggplot2)
library(scales)
library(dplyr)
library(readr)
library(Hmisc)
library(RColorBrewer)
library(grid)
library(gridExtra)


#Matrices - analysis

#load data
substitutions <- read.csv("Number_Of_Possibles.csv", header=FALSE)
edges <- read.csv("Average_Edges_Per_Node.csv", header=FALSE)
subs_edges <- merge(substitutions, edges, by.x="V1", by.y="V1")

Ideal <- read.csv("Parameters/Ideal_probs05.csv", header=FALSE)
Ideal$name <- "Equal Probabilities"

MutaT7 <- read.csv("Parameters/MutaT7_probs05.csv", header=FALSE)
MutaT7$name <- "MutaT7"

Mutazyme2 <- read.csv("Parameters/MutazymeII_probs05.csv", header=FALSE)
Mutazyme2$name <- "epPCR (Mutazyme II)"

OrthoRep <- read.csv("Parameters/OrthoRep_probs05.csv", header=FALSE)
OrthoRep$name <- "OrthoRep"

PACE <- read.csv("Parameters/PACE_probs05.csv", header=FALSE)
PACE$name <- "PACE"

TaqMn2 <- read.csv("Parameters/Taq_probs05.csv", header=FALSE)
TaqMn2$name <- "epPCR (Taq, Mn2+)"

TRIDENT <- read.csv("Parameters/TRIDENT_probs05.csv", header=FALSE)
TRIDENT$name <- "TRIDENT"

YeastIT1 <- read.csv("Parameters/YeastIT1_probs05.csv", header=FALSE)
YeastIT1$name <- "YeastIT1"

YeastIT2 <- read.csv("Parameters/YeastIT2_probs05.csv", header=FALSE)
YeastIT2$name <- "YeastIT2"

allprobs05 <- rbind(Ideal, MutaT7, Mutazyme2, OrthoRep, PACE, TaqMn2, TRIDENT, YeastIT1, YeastIT2)

Ideal <- read.csv("Parameters/Ideal_probsall.csv", header=FALSE)
Ideal$name <- "Equal Probabilities"

MutaT7 <- read.csv("Parameters/MutaT7_probsall.csv", header=FALSE)
MutaT7$name <- "MutaT7"

Mutazyme2 <- read.csv("Parameters/MutazymeII_probsall.csv", header=FALSE)
Mutazyme2$name <- "epPCR (Mutazyme II)"

OrthoRep <- read.csv("Parameters/OrthoRep_probsall.csv", header=FALSE)
OrthoRep$name <- "OrthoRep"

PACE <- read.csv("Parameters/PACE_probsall.csv", header=FALSE)
PACE$name <- "PACE"

TaqMn2 <- read.csv("Parameters/Taq_probsall.csv", header=FALSE)
TaqMn2$name <- "epPCR (Taq, Mn2+)"

TRIDENT <- read.csv("Parameters/TRIDENT_probsall.csv", header=FALSE)
TRIDENT$name <- "TRIDENT"

YeastIT1 <- read.csv("Parameters/YeastIT1_probsall.csv", header=FALSE)
YeastIT1$name <- "YeastIT1"

YeastIT2 <- read.csv("Parameters/YeastIT2_probsall.csv", header=FALSE)
YeastIT2$name <- "YeastIT2"

allprobsall <- rbind(Ideal, MutaT7, Mutazyme2, OrthoRep, PACE, TaqMn2, TRIDENT, YeastIT1, YeastIT2)


#plots
p1 <- ggplot(substitutions, aes(x=reorder(V1, -V2), y=V2, fill=c("#68AD7D")))
p1 <- p1 + geom_bar(stat="identity") + scale_fill_manual(values = c("#68AD7D")) + theme_grey(base_size=12)
p1 <- p1 + xlab("") + ylab("Substitutions with p≥0.05") 
p1 <- p1 + theme(legend.title = element_blank()) 
p1 <- p1 + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
#p1 <- p1 + scale_x_discrete(limits=c("Equal Probabilities", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "MutaT7", "PACE", "OrthoRep", "TRIDENT", "YeastIT1", "YeastIT2"))
#uncomment above for original order
p1 <- p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p1 <- p1 + theme(legend.position="none") 
p1 

p2 <- ggplot(edges, aes(x=reorder(V1, -V2), y=V2, fill=c("#68AD7D")))
p2 <- p2 + geom_bar(stat="identity") + scale_fill_manual(values = c("#68AD7D")) + theme_bw(base_size=12)
p2 <- p2 + xlab("") + ylab("Average number of edges/node") 
p2 <- p2 + theme(legend.title = element_blank()) 
p2 <- p2 + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
#p2 <- p2 + scale_x_discrete(limits=c("Equal Probabilities", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "MutaT7", "PACE", "OrthoRep", "TRIDENT", "YeastIT1", "YeastIT2"))
#uncomment above for original order
p2 <- p2 + theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
p2 <- p2 + theme(legend.position="none") 
p2

grid.arrange(arrangeGrob(p1, p2, widths=c(1,1), heights=c(1), nrow=1))

#graph with 2 axes for edges and substitutions
factor=1/21
p3 <- ggplot(subs_edges, aes(x = reorder(V1, -V2.x), y = V2.x)) +
  theme_bw(base_size=12)  +
  geom_bar(fill = "#68AD7D", stat = "identity") + geom_point(aes(x=V1, y=21*subs_edges$V2.y), colour="#DEE0E4", size=0) +
  scale_y_continuous(sec.axis = sec_axis(~.*factor, name="Average number of edges (p≥0.05)/node")) +
  #theme(
  #axis.title.y = element_text(color = c("#68AD7D"), size=12),
  #axis.title.y.right = element_text(color = c("#AE617B"), size=12)) + 
  xlab("") + ylab("Number of substitutions with p≥0.05") + theme(axis.text.x = element_text(angle = 50, hjust = 1)) 

p3

p4 <-ggplot(allprobs, aes(x=reorder(name, -(-V1)), y=V1)) +  theme_grey(base_size=12)  +
   geom_point(aes(color=allprobs$V1), position = position_jitter(width= 0.2, height = 0), size = 1) + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
   xlab("") + ylab("Probability of substitution") +
  scale_color_distiller(palette="YlGnBu", trans = "reverse") + theme(legend.title = element_blank()) + theme(legend.position="bottom") 
p4
grid.arrange(arrangeGrob(p1, p4, widths=c(1), heights=c(0.5,1), nrow=2))

p5 <-ggplot(allprobs05, aes(x=reorder(name, -(-V1)), y=V1)) + theme_grey(base_size=12) 
p5 <- p5 + geom_point(data=allprobs05, aes(x=reorder(name, -(-V1)), y=V1, color=allprobs05$V1), position = position_jitter(width= 0.2, height = 0), size = 1) + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") + ylab("Probability of substitution") +
  scale_color_distiller(palette="YlGnBu", trans = "reverse") + theme(legend.title = element_blank()) + theme(legend.position="bottom") 
p5 <- p5 + geom_boxplot(data=allprobsall, aes(x=reorder(name, -(-V1)), y=allprobsall$V1), outlier.shape=NA, alpha=0.1)
p5
grid.arrange(arrangeGrob(p1, p5, widths=c(1), heights=c(0.5,1), nrow=2))

#Directed graphs - analysis

graphs_stats <- read.csv("Codon_Graphs/Properties/stats.csv", header=TRUE)

graphs_stats$name[graphs_stats$file=="OrthoRep_Bias"] <- "OrthoRep"
graphs_stats$name[graphs_stats$file=="MutazymeII_Bias"] <- "epPCR (Mutazyme II)"
graphs_stats$name[graphs_stats$file=="Ideal_Bias"] <- "Equal Probabilities"
graphs_stats$name[graphs_stats$file=="PACE_Bias"] <- "PACE"
graphs_stats$name[graphs_stats$file=="MutaT7_Bias"] <- "(e)MutaT7"
graphs_stats$name[graphs_stats$file=="YeastIT2_Bias"] <- "YeastIT2"
graphs_stats$name[graphs_stats$file=="YeastIT1_Bias"] <- "YeastIT1"
graphs_stats$name[graphs_stats$file=="Taq_Bias"] <- "epPCR (Taq, Mn2+)"
graphs_stats$name[graphs_stats$file=="TRIDENT_Bias"] <- "TRIDENT"

#betweenness centrality - violin
p6 <-ggplot(graphs_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness)) + theme_grey(base_size=12) 
p6 <- p6 + geom_violin(data=graphs_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness), scale=3, draw_quantiles = TRUE) + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") + ylab("Betweenness centrality") +
  scale_color_distiller(palette="YlGnBu", trans = "reverse") + theme(legend.title = element_blank()) + theme(legend.position="bottom") 
#p6 <- p6 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness), outlier.shape=NA, alpha=0.1)
p6 <- p6 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
p6

#betweenness centrality
p7 <-ggplot(graphs_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness, color= betweenness)) + theme_grey(base_size=12) 
p7 <- p7 + geom_point(data=graphs_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness), color=c("#68AD7D"), position = position_jitter(width= 0.2, height = 0.007), size = 1)  + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") + ylab("Betweenness centrality") +
  scale_color_distiller(palette="Blues", trans = "reverse") + theme(legend.title = element_blank()) + theme(legend.position="bottom") 
p7 <- p7 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness), outlier.shape=NA, alpha=0.1)
p7 <- p7 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
p7 <- p7 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p7


#clustering coefficient
p8 <-ggplot(graphs_stats, aes(x=reorder(name, -(-clustering)), y=clustering)) + theme_grey(base_size=12) 
p8 <- p8 + geom_point(data=graphs_stats, aes(x=reorder(name, -(-clustering)), y=clustering),  color=c("#68AD7D"), position = position_jitter(width= 0.25, height = 0.01), size = 1)  + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") + ylab("Clustering coefficient") +
  scale_color_distiller(palette="Blues", trans = "reverse") + theme(legend.title = element_blank()) + theme(legend.position="none") 
p8 <- p8 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-clustering)), y=clustering), outlier.shape=NA, alpha=0.1)
p8 <- p8 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
p8 <- p8 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p8

#closeness centrality
p10 <-ggplot(graphs_stats, aes(x=reorder(name, -(-closeness_centrality)), y=closeness)) + theme_grey(base_size=12) 
p10 <- p10 + geom_point(data=graphs_stats, aes(x=reorder(name, -(-closeness)), y=closeness),  color=c("#68AD7D"), position = position_jitter(width= 0.25, height = 0.15), size = 1)  + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") + ylab("Closeness centrality") +
  theme(legend.title = element_blank()) + theme(legend.position="bottom") 
p10 <- p10 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-closeness)), y=closeness), outlier.shape=NA, alpha=0.1)
p10 <- p10 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
p10 <- p10 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p10

#degree centrality
p11 <-ggplot(graphs_stats, aes(x=reorder(name, -(-degree)), y=degree)) + theme_grey(base_size=12) 
p11 <- p11 + geom_point(data=graphs_stats, aes(x=reorder(name, -(-degree)), y=degree), color=c("#68AD7D"), position = position_jitter(width= 0.25, height = 0.008), size = 1)  
p11 <- p11 + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") + ylab("Degree centrality") +
  theme(legend.title = element_blank()) + theme(legend.position="bottom") 
p11 <- p11 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-degree)), y=degree), alpha=0.1, outlier.shape=NA)
p11 <- p11 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
p11 <- p11 + scale_color_distiller(palette="Blues", trans = "reverse") + theme(legend.title = element_blank()) + theme(legend.position="none") 

p11

#out-degree centrality
p12 <-ggplot(graphs_stats, aes(x=reorder(name, -(-out.degree)), y=out.egree)) + theme_grey(base_size=12) 
p12 <- p12 + geom_point(data=graphs_stats, aes(x=reorder(name, -(-out.degree)), y=out.degree), color=c("#68AD7D"), position = position_jitter(width= 0.25, height = 0.008), size = 1)  
p12 <- p12 + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") + ylab("Out-degree centrality") +
  theme(legend.title = element_blank()) + theme(legend.position="bottom") 
p12 <- p12 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-out.degree)), y=out.degree), alpha=0.1, outlier.shape=NA)
p12 <- p12 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
p12 <- p12 + scale_color_distiller(palette="Blues", trans = "reverse") + theme(legend.title = element_blank()) + theme(legend.position="none") 

p12


#higlighting data points

betweenness_highlight_stats <- graphs_stats %>% 
  filter(betweenness>=0.2 & name=="OrthoRep" | betweenness>=0.075 & name!="OrthoRep")
betweenness_rest_stats <- anti_join(graphs_stats, betweeness_highlight_stats)

closeness_highlight_stats <- graphs_stats %>% 
  filter(closeness>=9 & name=="epPCR (Mutazyme II)" | closeness>=6 & name=="epPCR (Taq, Mn2+)" | closeness>=2.4 & name=="YeastIT2")
closeness_rest_stats <- anti_join(graphs_stats, closeness_highlight_stats)

degree_highlight_stats <- graphs_stats %>% 
  filter(degree>=0.22 & name=="YeastIT2" | degree>=0.15 & name=="OrthoRep" | degree<0.07 & name=="TRIDENT" | degree<0.15 & name=="YeastIT2")
degree_rest_stats <- anti_join(graphs_stats, degree_highlight_stats)

outdegree_highlight_stats <- graphs_stats %>% 
  filter(out.degree>=0.125 & name=="YeastIT2" | out.degree>=0.075 & name=="OrthoRep" | out.degree<0.07 & name=="YeastIT2" | out.degree<0.01 & name=="(e)MutaT7")
outdegree_rest_stats <- anti_join(graphs_stats, outdegree_highlight_stats)



#betweenness centrality - highlighted

p7 <- ggplot(graphs_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness)) + theme_grey(base_size=12) 
p7 <- p7 + geom_point(data=betweenness_rest_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness), color=c("#68AD7D"), position = position_jitter(width= 0.2, height = 0.002), size = 1)  + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") + ylab("Betweenness centrality") + theme(legend.title = element_blank()) + theme(legend.position="none") 
p7 <- p7 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness), outlier.shape=NA, alpha=0.1)
p7 <- p7 + geom_point(data=betweenness_highlight_stats, aes(x=reorder(name, -(-betweenness)), y=betweenness, color=name), position = position_jitter(width= 0.2, height = 0.002), size = 1)
p7 <- p7 + scale_color_manual(name="legend", values=c("#AE617B", "#002255", "#002255"))
p7 <- p7 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
#p7 <- p7 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p7

#closeness centrality - highlighted

p10 <- ggplot(graphs_stats, aes(x=reorder(name, -(-closeness)), y=closeness)) + theme_grey(base_size=12) 
p10 <- p10 + geom_point(data=closeness_rest_stats, aes(x=reorder(name, -(-closeness)), y=closeness), color=c("#68AD7D"), position = position_jitter(width= 0.2, height = 0.1), size = 1)  + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") + ylab("Closeness centrality") + theme(legend.title = element_blank()) + theme(legend.position="none") 
p10 <- p10 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-closeness)), y=closeness), outlier.shape=NA, alpha=0.1)
p10 <- p10 + geom_point(data=closeness_highlight_stats, aes(x=reorder(name, -(-closeness)), y=closeness, color=name), position = position_jitter(width= 0.2, height = 0.1), size = 1)
p10 <- p10 + scale_color_manual(name="legend", values=c( "#002255", "#AE617B", "#AE617B"))
p10 <- p10 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
#p10 <- p10 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p10

#degree centrality - highlighted
p11 <-ggplot(graphs_stats, aes(x=reorder(name, -(-degree)), y=degree)) + theme_grey(base_size=12) 
p11 <- p11 + geom_point(data=degree_rest_stats, aes(x=reorder(name, -(-degree)), y=degree), color=c("#68AD7D"), position = position_jitter(width= 0.2, height = 0.004), size = 1)  
p11 <- p11 + theme(axis.text.x = element_text(angle = 50, hjust = 1)) + xlab("") + ylab("Degree centrality") +
  theme(legend.title = element_blank()) + theme(legend.position="none") 
p11 <- p11 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-degree)), y=degree), alpha=0.1, outlier.shape=NA)
p11 <- p11 + geom_point(data=degree_highlight_stats, aes(x=reorder(name, -(-degree)), y=degree, color=name), position = position_jitter(width= 0.2, height = 0.004), size = 1)
p11 <- p11 + scale_color_manual(name="legend", values=c("#AE617B", "#2c89a0","#A77EAD"))
p11 <- p11 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
p11 <- p11 + theme(legend.title = element_blank()) + theme(legend.position="none") 
p11

#out-degree centrality - highlighted
p12 <-ggplot(graphs_stats, aes(x=reorder(name, -(-out.degree)), y=out.degree)) + theme_grey(base_size=12) 
p12 <- p12 + geom_point(data=outdegree_rest_stats, aes(x=reorder(name, -(-out.degree)), y=out.degree), color=c("#68AD7D"), position = position_jitter(width= 0.2, height = 0.004), size = 1)  
p12 <- p12 + theme(axis.text.x = element_text(angle = 50, hjust = 1)) + xlab("") + ylab("Out-degree centrality") +
  theme(legend.title = element_blank()) + theme(legend.position="none") 
p12 <- p12 + geom_boxplot(data=graphs_stats, aes(x=reorder(name, -(-out.degree)), y=out.degree), alpha=0.1, outlier.shape=NA)
p12 <- p12 + geom_point(data=outdegree_highlight_stats, aes(x=reorder(name, -(-out.degree)), y=out.degree, color=name), position = position_jitter(width= 0.2, height = 0.004), size = 1)
p12 <- p12 + scale_color_manual(name="legend", values=c("#AE617B", "#002255","#AE617B"))
#need to change upper YeastIT2 group to dark blue manually
p12 <- p12 + scale_x_discrete(limits=c("Equal Probabilities", "PACE", "epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "YeastIT2", "OrthoRep", "YeastIT1", "TRIDENT", "(e)MutaT7"))
p12 <- p12 + theme(legend.title = element_blank()) + theme(legend.position="none") 
p12


grid.arrange(arrangeGrob(p7, p10, p12, widths=c(1,1,1), nrow=1))

grid.arrange(arrangeGrob(p1, p5, widths=c(1), heights=c(0.5,1), nrow=2))


