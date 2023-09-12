library(maftools)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library (ggpubr)
library(iNEXT)
library(actuar)
library(cowplot)
library(readr)

#load and process all individual maf files
maf_yeastit1 <- data.table::fread("allsamples_bc02.vcf.maf")
maf_yeastit2 <- data.table::fread("allsamples_bc03.vcf.maf")
yeastit1 = read.maf(maf=maf_yeastit1, vc_nonSyn = "Targeted_Region")
yeastit2 = read.maf(maf=maf_yeastit2, vc_nonSyn = "Targeted_Region")
yeastit1
yeastit2

maf_yeastit1_cds <- data.table::fread("cdsonly_allsamples_bc02_head.vcf.maf")
yeastit1_cds = read.maf(maf=maf_yeastit1_cds, vc_nonSyn = "Targeted_Region")

maf_yeastit2_cds <- data.table::fread("cdsonly_allsamples_bc03_head.vcf.maf")
yeastit2_cds = read.maf(maf=maf_yeastit2_cds, vc_nonSyn = "Targeted_Region")


getSampleSummary(yeastit1_cds)
getSampleSummary(yeastit2_cds)
yeastit1_cds_mutfreq <- getSampleSummary(yeastit1_cds)
yeastit2_cds_mutfreq <- getSampleSummary(yeastit2_cds)

maf_yeastit1_cds_unique <- data.table::fread("cdsonly_allsamples_bc02_unique_head.vcf.maf")
yeastit1_cds_unique = read.maf(maf=maf_yeastit1_cds_unique, vc_nonSyn = "Targeted_Region")
maf_yeastit2_cds_unique <- data.table::fread("cdsonly_allsamples_bc03_unique_head.vcf.maf")
yeastit2_cds_unique = read.maf(maf=maf_yeastit2_cds_unique, vc_nonSyn = "Targeted_Region")

yeastit1_cds_unique_mutfreq <- getSampleSummary(yeastit1_cds_unique)
yeastit2_cds_unique_mutfreq <- getSampleSummary(yeastit2_cds_unique)

yeastit1_cds_titv_unique = titv(maf = yeastit1_cds_unique, plot = TRUE, useSyn = TRUE, file=NULL)
yeastit2_cds_titv_unique = titv(maf = yeastit2_cds_unique, plot = TRUE, useSyn = TRUE, file=NULL)


#mutation spectra analysis
mutation_types <- c('T>C','C>T','T>G', 'T>A', 'C>A', 'C>G')

#mutation spectra - cds only, all UMI clusters
yeastit1_cds_titv = titv(maf = yeastit1_cds, plot = TRUE, useSyn = TRUE, file=NULL)
yeastit2_cds_titv = titv(maf = yeastit2_cds, plot = TRUE, useSyn = TRUE, file=NULL)
mutation_values_yit1_cds <- c(sum(yeastit1_cds_titv$raw.counts$`T>C`),  sum(yeastit1_cds_titv$raw.counts$`C>T`),  sum(yeastit1_cds_titv$raw.counts$`T>G`),  sum(yeastit1_cds_titv$raw.counts$`T>A`),  sum(yeastit1_cds_titv$raw.counts$`C>A`),  sum(yeastit1_cds_titv$raw.counts$`C>G`))
mutation_values_yit2_cds <- c(sum(yeastit2_cds_titv$raw.counts$`T>C`),  sum(yeastit2_cds_titv$raw.counts$`C>T`),  sum(yeastit2_cds_titv$raw.counts$`T>G`),  sum(yeastit2_cds_titv$raw.counts$`T>A`),  sum(yeastit2_cds_titv$raw.counts$`C>A`),  sum(yeastit2_cds_titv$raw.counts$`C>G`))
mutation_spectrum_cds <- data.frame(mutation_types, mutation_values_yit1_cds, mutation_values_yit2_cds)  

#mutation spectra analysis - full-length, all UMI clusters
yeastit1_titv = titv(maf = yeastit1, plot = TRUE, useSyn = TRUE, file=NULL)
yeastit2_titv = titv(maf = yeastit2, plot = TRUE, useSyn = TRUE, file=NULL)
mutation_values_yit1 <- c(sum(yeastit1_titv$raw.counts$`T>C`),  sum(yeastit1_titv$raw.counts$`C>T`),  sum(yeastit1_titv$raw.counts$`T>G`),  sum(yeastit1_titv$raw.counts$`T>A`),  sum(yeastit1_titv$raw.counts$`C>A`),  sum(yeastit1_titv$raw.counts$`C>G`))
mutation_values_yit2 <- c(sum(yeastit2_titv$raw.counts$`T>C`),  sum(yeastit2_titv$raw.counts$`C>T`),  sum(yeastit2_titv$raw.counts$`T>G`),  sum(yeastit2_titv$raw.counts$`T>A`),  sum(yeastit2_titv$raw.counts$`C>A`),  sum(yeastit2_titv$raw.counts$`C>G`))
mutation_spectrum <- data.frame(mutation_types, mutation_values_yit1, mutation_values_yit2)  

#mutation spectra analysis - cds only, unique cluster analysis
mutation_values_yit1_cds_unique <- c(sum(yeastit1_cds_titv_unique$raw.counts$`T>C`),  sum(yeastit1_cds_titv_unique$raw.counts$`C>T`),  sum(yeastit1_cds_titv_unique$raw.counts$`T>G`),  sum(yeastit1_cds_titv_unique$raw.counts$`T>A`),  sum(yeastit1_cds_titv_unique$raw.counts$`C>A`),  sum(yeastit1_cds_titv_unique$raw.counts$`C>G`))
mutation_values_yit2_cds_unique <- c(sum(yeastit2_cds_titv_unique$raw.counts$`T>C`),  sum(yeastit2_cds_titv_unique$raw.counts$`C>T`),  sum(yeastit2_cds_titv_unique$raw.counts$`T>G`),  sum(yeastit2_cds_titv_unique$raw.counts$`T>A`),  sum(yeastit2_cds_titv_unique$raw.counts$`C>A`),  sum(yeastit2_cds_titv_unique$raw.counts$`C>G`))
mutation_spectrum_cds_unique <- data.frame(mutation_types, mutation_values_yit1_cds_unique, mutation_values_yit2_cds_unique)  


#mutation rate analysis - cds only 

#add lines for clones that had 0 mutation to files with non-unique clusters:

for (i in 7199:10503) {
  yeastit1_cds_mutfreq <- rbind(yeastit1_cds_mutfreq, data.frame(Targeted_Region = 0), fill=TRUE)
}
for (i in 881:4336) {
  yeastit2_cds_mutfreq <- rbind(yeastit2_cds_mutfreq, data.frame(Targeted_Region = 0), fill=TRUE)
}

#compare average mutation no. in all UMI clusters (with >= 1 point mutation) vs. unique clusters 
mean(yeastit1_cds_mutfreq$Targeted_Region[0:7198])
mean(yeastit1_cds_unique_mutfreq$Targeted_Region)
mean(yeastit2_cds_mutfreq$Targeted_Region[0:880])
mean(yeastit2_cds_unique_mutfreq$Targeted_Region)

#as the values above are similar, average mutation frequencies in yeastit1_cds_mutfreq 
#(including those with 0 point mutations) can be considered a good approximation for overall mutation frequency

mean(yeastit1_cds_mutfreq$Targeted_Region) #mutation rate - yeastIT1 /gene/3day
mean(yeastit2_cds_mutfreq$Targeted_Region) # mutation rate - yeastIT2 /gene/3day


#analysis - variant vs. cluster count (full-length)

yeast1_clusterlist <- data.frame(maf_yeastit1$Start_Position, maf_yeastit1$Tumor_Seq_Allele2, maf_yeastit1$Tumor_Sample_Barcode)
yeast1_clusterlist <- aggregate(data=yeast1_clusterlist, cbind(maf_yeastit1.Start_Position, maf_yeastit1.Tumor_Seq_Allele2)~maf_yeastit1.Tumor_Sample_Barcode, FUN=paste)
colnames(yeast1_clusterlist) = c("ClusterID","SNP position", "Mutation")
yeast1_clusterlist2 <- rename(count(yeast1_clusterlist, `SNP position`, `Mutation`), Frequency = n)
yeast1_clusterlist2_ordered <- yeast1_clusterlist2[order(yeast1_clusterlist2$Frequency),]
yeast1_clusterlist2_ordered$ClusterID <- c(1:669)
yeast1_clusterlist2_ordered$mutno <- sapply(yeast1_clusterlist2_ordered$`SNP position`, length)


yeast2_clusterlist <- data.frame(maf_yeastit2$Start_Position, maf_yeastit2$Tumor_Seq_Allele2, maf_yeastit2$Tumor_Sample_Barcode)
yeast2_clusterlist <- aggregate(data=yeast2_clusterlist, cbind(maf_yeastit2.Start_Position, maf_yeastit2.Tumor_Seq_Allele2)~maf_yeastit2.Tumor_Sample_Barcode, FUN=paste)
colnames(yeast2_clusterlist) = c("ClusterID","SNP position", "Mutation")
yeast2_clusterlist2 <- rename(count(yeast2_clusterlist, `SNP position`, `Mutation`), Frequency = n)
yeast2_clusterlist2_ordered <- yeast2_clusterlist2[order(yeast2_clusterlist2$Frequency),]
yeast2_clusterlist2_ordered$ClusterID <- c(1:309)
yeast2_clusterlist2_ordered$mutno <- sapply(yeast2_clusterlist2_ordered$`SNP position`, length)


#population diversity estimation
yeast1_clusterlist2_ordered_add0 <- rbind(yeast1_clusterlist2_ordered, c(0,0,3305,670,0))
yeast2_clusterlist2_ordered_add0 <- rbind(yeast2_clusterlist2_ordered, c(0,0,3556,310,0))

mutno=1/60
p11 <- ggplot(yeast1_clusterlist2_ordered_add0, aes(x = ClusterID, y = Frequency)) +
  theme_bw(base_size=12)  +
  geom_bar(fill = "#68AD7D", stat = "identity") + geom_point(aes(x=ClusterID, y=60*yeast1_clusterlist2_ordered_add0$mutno), colour="#DEE0E4", size=1) +
  scale_y_continuous(sec.axis = sec_axis(~.*mutno, name="Mutation count")) +
  theme(
    axis.title.y = element_text(color = c("#68AD7D"), size=12),
    axis.title.y.right = element_text(color = c("#8C838D"), size=12)
  ) + xlab("Variant ID") + ylab("Cluster count")

p11

p12 <- ggplot(yeast2_clusterlist2_ordered_add0, aes(x = ClusterID, y = Frequency)) +
  theme_bw(base_size=12) +
  geom_bar(fill = "#68AD7D", stat = "identity") + geom_point(aes(x=ClusterID, y=60*mutno), colour="#DEE0E4", size=1) +
  scale_y_continuous(sec.axis = sec_axis(~.*mutno, name="Mutation count")) +
  theme(
    axis.title.y = element_text(color = c("#68AD7D"), size=12),
    axis.title.y.right = element_text(color = c("#8C838D"), size=12))  +
  xlab("Variant ID") + ylab("Cluster count")

p12

grid.arrange(arrangeGrob(p11, p12, widths=c(1,1), nrow=1)) 


yeast1_clusterlist2_ordered_freqtable <- subset(yeast1_clusterlist2_ordered_add0, select=c("Frequency"))
yeast2_clusterlist2_ordered_freqtable <- subset(yeast2_clusterlist2_ordered_add0, select=c("Frequency"))

yeast1_iNEXT <- iNEXT(yeast1_clusterlist2_ordered_freqtable, q=c(0), datatype="abundance", endpoint=200000, knots =30)
yeast2_iNEXT <- iNEXT(yeast2_clusterlist2_ordered_freqtable, q=c(0), datatype="abundance", endpoint=100000, knots=30)

#Coverage-based R/E sampling curves: iNEXT computes diversity estimates for
# rarefied and extrapolated samples with sample completeness (as measured by
# sample coverage) up to the coverage value of double the reference sample size (bydefault) or 
# a user-specified coverage. This type of sampling curve plots the
# diversity estimates with respect to sample coverage. 
#yeast1_iNEXT2 <- estimateD(yeast1_clusterlist2_ordered_freqtable, q = c(0,1,2), datatype = "abundance", base="coverage", level=0.99, conf=0.95)

yeast1_iNEXT2 <- estimateD(yeast1_clusterlist2_ordered_freqtable, q = c(0,1,2), datatype = "abundance", base="size", level=200000)
yeast2_iNEXT2 <- estimateD(yeast2_clusterlist2_ordered_freqtable, q = c(0,1,2), datatype = "abundance", base="size", level=200000)



#p12 and p13 show estimation of library difersity after 1st transformation to E.coli where approximately 10k 
#colonies were obtained for each variant

p13 <- ggiNEXT(yeast1_iNEXT, type=1, se=TRUE) + theme_bw(base_size=12) + theme(legend.position="none")
p13 <- p13 + xlab("Sample size")
p13 <- p13 + ylab("Estimated number of unique variants")
p13


p14 <- ggiNEXT(yeast2_iNEXT, type=1, se=TRUE) + theme_bw(base_size=12) + theme(legend.position="none")
p14 <- p14 + xlab("Sample size")
p14 <- p14 + ylab("Estimated number of unique variants")
p14

grid.arrange(arrangeGrob(p13, p14, widths=c(1,1), nrow=1)) 



#mutation spectra plots
#boxplots - mutation spectra - unique variants, cds only

#yit1-format data for boxplots
yit1_cds_tc <- data.frame("AT>GC", yeastit1_cds_titv_unique$raw.counts$`T>C`)
colnames(yit1_cds_tc) <- c("mutation", "events")

yit1_cds_ct <- data.frame("GC>AT", yeastit1_cds_titv_unique$raw.counts$`C>T`)
colnames(yit1_cds_ct) <- c("mutation", "events")

yit1_cds_tg <- data.frame("AT>CG", yeastit1_cds_titv_unique$raw.counts$`T>G`)
colnames(yit1_cds_tg) <- c("mutation", "events")

yit1_cds_ta <- data.frame("AT>TA", yeastit1_cds_titv_unique$raw.counts$`T>A`)
colnames(yit1_cds_ta) <- c("mutation", "events")

yit1_cds_ca <- data.frame("GC>TA", yeastit1_cds_titv_unique$raw.counts$`C>A`)
colnames(yit1_cds_ca) <- c("mutation", "events")

yit1_cds_cg <- data.frame("GC>CG", yeastit1_cds_titv_unique$raw.counts$`C>G`)
colnames(yit1_cds_cg) <- c("mutation", "events")

mutation_list_yit1_cds <- rbind(yit1_cds_tc, yit1_cds_ct, yit1_cds_tg, yit1_cds_ta, yit1_cds_ca, yit1_cds_cg)

mutation_list_yit1_cds2 <-
  mutation_list_yit1_cds %>%
  group_by(mutation) %>%
  mutate(outlier = events > median(events) + IQR(events) * 1.5) %>%
  ungroup

yit1_cds_ti <- data.frame("Ti", (yeastit1_cds_titv_unique$raw.counts$`C>T`+yeastit1_cds_titv_unique$raw.counts$`T>C`))
colnames(yit1_cds_ti) <- c("mutation", "events")

yit1_cds_tv <- data.frame("Tv", (yeastit1_cds_titv_unique$raw.counts$`C>G`+yeastit1_cds_titv_unique$raw.counts$`C>A`+yeastit1_cds_titv_unique$raw.counts$`T>A`+yeastit1_cds_titv_unique$raw.counts$`T>G`))
colnames(yit1_cds_tv) <- c("mutation", "events")

yit1_cds_all <- data.frame("All", (yeastit1_cds_titv_unique$raw.counts$`C>T`+yeastit1_cds_titv_unique$raw.counts$`T>C`+yeastit1_cds_titv_unique$raw.counts$`C>G`+yeastit1_cds_titv_unique$raw.counts$`C>A`+yeastit1_cds_titv_unique$raw.counts$`T>A`+yeastit1_cds_titv_unique$raw.counts$`T>G`))
colnames(yit1_cds_all) <- c("mutation", "events")

mutation_list_yit1_cds_titv <- rbind(yit1_cds_ti,yit1_cds_tv,yit1_cds_all)

mutation_list_yit1_cds_titv2 <-
  mutation_list_yit1_cds_titv %>%
  group_by(mutation) %>%
  mutate(outlier = events > median(events) + IQR(events) * 1.5) %>%
  ungroup

#yit2-format data for boxplots
yit2_cds_tc <- data.frame("AT>GC", yeastit2_cds_titv_unique$raw.counts$`T>C`)
colnames(yit2_cds_tc) <- c("mutation", "events")

yit2_cds_ct <- data.frame("GC>AT", yeastit2_cds_titv_unique$raw.counts$`C>T`)
colnames(yit2_cds_ct) <- c("mutation", "events")

yit2_cds_tg <- data.frame("AT>CG", yeastit2_cds_titv_unique$raw.counts$`T>G`)
colnames(yit2_cds_tg) <- c("mutation", "events")

yit2_cds_ta <- data.frame("AT>TA", yeastit2_cds_titv_unique$raw.counts$`T>A`)
colnames(yit2_cds_ta) <- c("mutation", "events")

yit2_cds_ca <- data.frame("GC>TA", yeastit2_cds_titv_unique$raw.counts$`C>A`)
colnames(yit2_cds_ca) <- c("mutation", "events")

yit2_cds_cg <- data.frame("GC>CG", yeastit2_cds_titv_unique$raw.counts$`C>G`)
colnames(yit2_cds_cg) <- c("mutation", "events")


mutation_list_yit2_cds <- rbind(yit2_cds_tc, yit2_cds_ct, yit2_cds_tg, yit2_cds_ta, yit2_cds_ca, yit2_cds_cg)

mutation_list_yit2_cds2 <-
  mutation_list_yit2_cds %>%
  group_by(mutation) %>%
  mutate(outlier = events > median(events) + IQR(events) * 1.5) %>%
  ungroup


yit2_cds_ti <- data.frame("Ti", (yeastit2_cds_titv_unique$raw.counts$`C>T`+yeastit2_cds_titv_unique$raw.counts$`T>C`))
colnames(yit2_cds_ti) <- c("mutation", "events")

yit2_cds_tv <- data.frame("Tv", (yeastit2_cds_titv_unique$raw.counts$`C>G`+yeastit2_cds_titv_unique$raw.counts$`C>A`+yeastit2_cds_titv_unique$raw.counts$`T>A`+yeastit2_cds_titv_unique$raw.counts$`T>G`))
colnames(yit2_cds_tv) <- c("mutation", "events")

yit2_cds_all <- data.frame("All", (yeastit2_cds_titv_unique$raw.counts$`C>T`+yeastit2_cds_titv_unique$raw.counts$`T>C`+yeastit2_cds_titv_unique$raw.counts$`C>G`+yeastit2_cds_titv_unique$raw.counts$`C>A`+yeastit2_cds_titv_unique$raw.counts$`T>A`+yeastit2_cds_titv_unique$raw.counts$`T>G`))
colnames(yit2_cds_all) <- c("mutation", "events")

mutation_list_yit2_cds_titv <- rbind(yit2_cds_ti,yit2_cds_tv,yit2_cds_all)

mutation_list_yit2_cds_titv2 <-
  mutation_list_yit2_cds_titv %>%
  group_by(mutation) %>%
  mutate(outlier = events > median(events) + IQR(events) * 1.5) %>%
  ungroup

#plots box plots -  mutation counts by type

#yit1
p3 <- ggplot(mutation_list_yit1_cds2, aes(fill=factor(`mutation`, levels=c("AT>GC", "GC>AT", "AT>CG", "AT>TA", "GC>TA", "GC>CG")), x=mutation, y=events)) + geom_boxplot(outlier.shape=NA)
p3 <- p3 + xlab("Mutation type") + ylab("Count")
p3 <- p3 + scale_x_discrete(limits=c("AT>GC", "GC>AT", "AT>CG", "AT>TA", "GC>TA", "GC>CG"))
p3 <- p3 + theme_bw(base_size=11) + theme(axis.text.x = element_text(angle = 50, hjust = 1))
p3 <- p3 + scale_fill_manual(values = c("#68AD7D", "#97D2A9", "#A77EAD", "#7c6f91ff", "#8C838D",  "#AE617B"))
p3 <- p3 + geom_point(aes(colour=factor(mutation, levels=c("AT>GC", "GC>AT", "AT>CG", "AT>TA", "GC>TA", "GC>CG"))), data=function(x) dplyr::filter_(x, ~ outlier), position = 'jitter', alpha=0.7) 
p3 <- p3 + scale_colour_manual(values = c("#68AD7D", "#97D2A9", "#A77EAD", "#7c6f91ff", "#8C838D",  "#AE617B"))
p3 <- p3 + theme(legend.title = element_blank()) + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
p3 <- p3 + theme(legend.position = "none")
p3


p4 <- ggplot(mutation_list_yit1_cds_titv2, aes(fill=factor(`mutation`, levels=c("Ti", "Tv", "All")), x=mutation, y=events)) + geom_boxplot(outlier.shape=NA)
p4 <- p4 + xlab("Mutation type") + ylab("Count")
p4 <- p4 + scale_x_discrete(limits=c("Ti", "Tv", "All"))
p4 <- p4 + theme_bw(base_size=11) 
p4 <- p4 + scale_fill_manual(values = c("#68AD7D", "#A77EAD", "#DEE0E4"))
p4 <- p4 + geom_point(aes(colour=factor(mutation,  levels=c("Ti", "Tv", "All"))), data=function(x) dplyr::filter_(x, ~ outlier), position = 'jitter', alpha=0.7) 
p4 <- p4 + scale_colour_manual(values = c("#68AD7D", "#A77EAD", "#DEE0E4"))
p4 <- p4 + theme(legend.title = element_blank()) + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
p4 <- p4 + theme(legend.position = "none")
p4

#yit2

p5 <- ggplot(mutation_list_yit2_cds2, aes(fill=factor(`mutation`, levels=c("AT>GC", "GC>AT", "AT>CG", "AT>TA", "GC>TA", "GC>CG")), x=mutation, y=events)) + geom_boxplot(outlier.shape=NA)
p5 <- p5 + xlab("Mutation type") + ylab("Count")
p5 <- p5 + scale_x_discrete(limits=c("AT>GC", "GC>AT", "AT>CG", "AT>TA", "GC>TA", "GC>CG"))
p5 <- p5 + theme_bw(base_size=11) + theme(axis.text.x = element_text(angle = 50, hjust = 1))
p5 <- p5 + scale_fill_manual(values = c("#68AD7D", "#97D2A9", "#A77EAD", "#7c6f91ff", "#8C838D",  "#AE617B"))
p5 <- p5 + geom_point(aes(colour=factor(mutation, levels=c("AT>GC", "GC>AT", "AT>CG", "AT>TA", "GC>TA", "GC>CG"))), data=function(x) dplyr::filter_(x, ~ outlier), position = 'jitter', alpha=0.7) 
p5 <- p5 + scale_colour_manual(values = c("#68AD7D", "#97D2A9", "#A77EAD", "#7c6f91ff", "#8C838D",  "#AE617B"))
p5 <- p5 + theme(legend.title = element_blank()) + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
p5 <- p5 + theme(legend.position = "none")
p5


p6 <- ggplot(mutation_list_yit2_cds_titv2, aes(fill=factor(`mutation`, levels=c("Ti", "Tv", "All")), x=mutation, y=events)) + geom_boxplot(outlier.shape=NA)
p6 <- p6 + xlab("Mutation type") + ylab("Count")
p6 <- p6 + scale_x_discrete(limits=c("Ti", "Tv", "All"))
p6 <- p6 + theme_bw(base_size=11) 
p6 <- p6 + geom_point(aes(colour=factor(mutation, levels=c("Ti", "Tv", "All"))), data=function(x) dplyr::filter_(x, ~ outlier), position = 'jitter', alpha=0.7) 
p6 <- p6 + scale_colour_manual(values = c("#68AD7D", "#A77EAD", "#DEE0E4"))
p6 <- p6 + theme(legend.title = element_blank()) + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
p6 <- p6 + theme(legend.position = "none")
p6


#plot histograms: Frequency distribution of UMI-linked sequence variants with different number of mutations
p1 <- ggplot(yeastit1_cds_mutfreq, aes(x=yeastit1_cds_mutfreq$Targeted_Region)) + geom_histogram(binwidth=2, fill=c("#DEE0E4")) + theme_bw(base_size=11)
p1 <- p1 + xlab("Number of mutations") + ylab("Count")
p1 <- p1 + geom_vline(xintercept=mean(yeastit1_cds_mutfreq$Targeted_Region), lwd=1, linetype=2, color="black")
p1

p2 <- ggplot(yeastit2_cds_mutfreq, aes(x=yeastit2_cds_mutfreq$Targeted_Region)) + geom_histogram(binwidth=1, fill=c("#DEE0E4")) + theme_bw(base_size=11)
p2 <- p2 + xlab("Number of mutations") + ylab("Count")
p2 <- p2 + geom_vline(xintercept=mean(yeastit2_cds_mutfreq$Targeted_Region), lwd=1, linetype=2, color="black")

p2

grid.arrange(arrangeGrob(p1, p2, p4, p6, p3, p5, widths=c(1,1), heights=c(1,1,1.3), nrow=3))

plot_grid(p1, p2, p4, p6, p3, p5, align = "v", nrow = 3, rel_heights = c(1, 1, 1.2))


#plot histograms: Frequency distribution of unique sequence variants with 
#different number of mutations (excl.0)
#QC only - to confirm that distributions align well

p8 <- ggplot(yeastit1_cds_unique_mutfreq, aes(x=yeastit1_cds_unique_mutfreq$Targeted_Region)) + geom_histogram(binwidth=1, fill=c("#DEE0E4")) + theme_bw(base_size=12)
p8 <- p8 + xlab("Number of mutations") + ylab("Count")
p8

p9 <- ggplot(yeastit1_cds_mutfreq, aes(x=yeastit1_cds_mutfreq$Targeted_Region)) + geom_histogram(binwidth=2, fill=c("#DEE0E4")) + theme_bw(base_size=12)
p9 <- p9 + xlab("Number of mutations") + ylab("Count")
p9

p10 <- ggplot(yeastit2_cds_unique_mutfreq, aes(x=yeastit2_cds_unique_mutfreq$Targeted_Region)) + geom_histogram(binwidth=1, fill=c("#DEE0E4")) + theme_bw(base_size=12)
p10 <- p10 + xlab("Number of mutations") + ylab("Count") + scale_x_continuous(breaks=c(0,1,2,3,4,5))
p10

p11 <- ggplot(yeastit2_cds_mutfreq, aes(x=yeastit2_cds_mutfreq$Targeted_Region)) + geom_histogram(binwidth=1, fill=c("#DEE0E4")) + theme_bw(base_size=12)
p11 <- p11 + xlab("Number of mutations") + ylab("Count") + scale_x_continuous(breaks=c(0,1,2,3,4,5))
p11


grid.arrange(arrangeGrob(p9, p8, p11, p10, widths=c(1,1), heights=c(1,1), nrow=2)) #comparison of mutational frequency distribution, YIT1 vs YIT2, all UMIs vs unique variants


#generate frequency SNP tables for lollipop plots

maf_yeastit1_unique <- data.table::fread("allsamples_bc02_unique_head.vcf.maf")
maf_yeastit2_unique <- data.table::fread("allsamples_bc03_unique_head.vcf.maf")

yeast1_allSNPlist_unique <- data.frame(maf_yeastit1_unique$Start_Position, maf_yeastit1_unique$Tumor_Seq_Allele1, maf_yeastit1_unique$Tumor_Seq_Allele2)
yeast2_allSNPlist_unique <- data.frame(maf_yeastit2_unique$Start_Position, maf_yeastit2_unique$Tumor_Seq_Allele1, maf_yeastit2_unique$Tumor_Seq_Allele2)

yeast1_allSNPlist_unique <- yeast1_allSNPlist_unique[!(yeast1_allSNPlist_unique$maf_yeastit1_unique.Tumor_Seq_Allele2=="-"),] #remove non SNP mutations
yeast1_allSNPlist_unique <- yeast1_allSNPlist_unique[!(yeast1_allSNPlist_unique$maf_yeastit1_unique.Tumor_Seq_Allele1=="-"),] #remove non SNP mutations
yeast1_allSNPlist_unique$`Mutation` <- paste(yeast1_allSNPlist_unique$maf_yeastit1_unique.Tumor_Seq_Allele1, ">", yeast1_allSNPlist_unique$maf_yeastit1_unique.Tumor_Seq_Allele2)
yeastit1_allSNP_unique_counts <- rename(count(yeast1_allSNPlist_unique, maf_yeastit1_unique.Start_Position, `Mutation`), Frequency = n)
yeastit1_allSNP_unique_counts$`Nt Position` <- c(yeastit1_allSNP_unique_counts$maf_yeastit1_unique.Start_Position - 251) #shift position numbers so that 1st nucleotide of start codon = 0


yeast2_allSNPlist_unique <- yeast2_allSNPlist_unique[!(yeast2_allSNPlist_unique$maf_yeastit2_unique.Tumor_Seq_Allele2=="-"),] #remove non SNP mutations
yeast2_allSNPlist_unique <- yeast2_allSNPlist_unique[!(yeast2_allSNPlist_unique$maf_yeastit2_unique.Tumor_Seq_Allele1=="-"),] #remove non SNP mutations
yeast2_allSNPlist_unique$`Mutation` <- paste(yeast2_allSNPlist_unique$maf_yeastit2_unique.Tumor_Seq_Allele1, ">", yeast2_allSNPlist_unique$maf_yeastit2_unique.Tumor_Seq_Allele2)
yeastit2_allSNP_unique_counts <- rename(count(yeast2_allSNPlist_unique, maf_yeastit2_unique.Start_Position, `Mutation`), Frequency = n)
yeastit2_allSNP_unique_counts$`Nt Position` <- c(yeastit2_allSNP_unique_counts$maf_yeastit2_unique.Start_Position - 251) #shift position numbers so that 1st nucleotide of start codon = 0

yeastit1_allSNP_unique_counts$class <- "YeastIT1"
colnames(yeastit1_allSNP_unique_counts)=c("Start_Position", "Mutation", "Frequency", "Nt position", "class")
yeastit2_allSNP_unique_counts$class <- "YeastIT2"
colnames(yeastit2_allSNP_unique_counts)=c("Start_Position", "Mutation", "Frequency", "Nt position", "class")
merged_allSNP_unique_counts <- rbind(yeastit1_allSNP_unique_counts, yeastit2_allSNP_unique_counts)
merged_allSNP_unique_counts <- merged_allSNP_unique_counts[(merged_allSNP_unique_counts$`Mutation`=="A > G" | merged_allSNP_unique_counts$`Mutation`=="T > C" | merged_allSNP_unique_counts$`Mutation`=="G > A" | merged_allSNP_unique_counts$`Mutation`=="C > T" | merged_allSNP_unique_counts$`Mutation`=="A > C" | merged_allSNP_unique_counts$`Mutation`=="T > G" | merged_allSNP_unique_counts$`Mutation`=="A > T" | merged_allSNP_unique_counts$`Mutation`=="T > A" | merged_allSNP_unique_counts$`Mutation`=="G > T" | merged_allSNP_unique_counts$`Mutation`=="C > A" | merged_allSNP_unique_counts$`Mutation`=="G > C" | merged_allSNP_unique_counts$`Mutation`=="C > G"),] #clean up mutations further - removes all but SNP mutations, should not be used if many DNP mutations present

#plot lollipop plots

p7 <- ggplot(merged_allSNP_unique_counts, aes(x=`Nt position`, y=Frequency)) +
  facet_wrap(~class, ncol=1) +
  geom_segment( aes(x=`Nt position`, xend=`Nt position`, y=0.5, yend=Frequency), color="#DEE0E4") +
  geom_point(aes(colour=factor(Mutation, levels=c("A > G", "T > C", "G > A", "C > T", "A > C", "T > G", "A > T", "T > A", "G > T", "C > A", "G > C", "C > G")), shape=factor(Mutation, levels=c("A > G", "T > C", "G > A", "C > T", "A > C", "T > G", "A > T", "T > A", "G > T", "C > A", "G > C", "C > G"))), size=2) +
  theme_bw(base_size=12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Position relative to start codon") +
  ylab("Mutation frequency") + scale_y_log10()
p7 <- p7 + scale_colour_manual(values = c("#68AD7D","#68AD7D", "#97D2A9","#97D2A9", "#A77EAD","#A77EAD", "#7c6f91ff", "#7c6f91ff",  "#8C838D", "#8C838D", "#AE617B", "#AE617B"))
p7 <- p7 + scale_shape_manual(values=c(16,17,17,16,16,17,16,17,17,16,17,16))
p7 <- p7 + theme(legend.title = element_blank()) + theme(plot.background = element_rect(fill = "transparent"),  legend.box.background = element_rect(fill = "transparent", colour=NA))
p7 <- p7 + theme(legend.position="top") + guides(colour=guide_legend(nrow =1))
p7


#mutation spectrum analysis of all methods


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
  scale_x_discrete(limits=c("epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "EvolVR", "MutaT7", "eMutaT7", "PACE", "OrthoRep", "TRIDENT", "YeastIT1", "YeastIT2"))
p17 <- p17 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
#p17 <- p17 + theme(legend.position="top") + guides(colour=guide_legend(nrow =1))
p17 <- p17 + theme(legend.position="none") 
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
p18 <- p18 + scale_x_discrete(limits=c("epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "EvolVR", "MutaT7", "eMutaT7", "PACE", "OrthoRep", "TRIDENT", "YeastIT", "YeastIT2"))
p18 <- p18 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
#p18 <- p18 + theme(legend.position="top") + guides(colour=guide_legend(nrow =1))
p18 <- p18 + theme(legend.position="none") 
p18


#rates$`Mutation rate` <- as.factor(rates$`Mutation rate`)
p19 <- ggplot(rates, aes(x=Method, y=rates$`Mutation rate`, fill=c("#68AD7D")))
p19 <- p19 + geom_bar(stat="identity") + scale_fill_manual(values = c("#68AD7D")) + theme_bw(base_size=12)
p19 <- p19 + xlab("") + ylab(expression("Mutation rate"~(kb^{-1}~day^{-1}))) 
p19 <- p19 + theme(legend.title = element_blank()) 
p19 <- p19 + theme(plot.background = element_rect(fill = "transparent",colour = NA),  legend.box.background = element_rect(fill = "transparent", colour=NA))
p19 <- p19 + scale_x_discrete(limits=c("epPCR (Mutazyme II)", "epPCR (Taq, Mn2+)", "EvolVR", "MutaT7", "eMutaT7", "PACE", "OrthoRep", "TRIDENT", "YeastIT1", "YeastIT2"))
p19 <- p19 + theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
p19 <- p19 + theme(legend.position="none") 
p19   

p20 <- ggarrange(p17,p18,p19,nrow=3,ncol=1, heights=c(0.8,0.9,1.05))
p20

plot_grid(p1, p2, p17, p4, p6, p18, p3, p5, p19, align = "v", nrow = 3, rel_heights = c(1, 1, 1.2), rel_widths = c(1, 1, 2))

