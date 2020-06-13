setwd("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/1_ATACseq/2_output_abundance")

# READING BED FILES OF G BANDS

#############################
###         Gneg1         ###
#############################

ATACSeq_Gneg1 <- read.delim("ATACSeq_Gneg1.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gneg1 <- ATACSeq_Gneg1[,c(-4,-5,-6)]
colnames(ATACSeq_Gneg1) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gneg1$Abundance <- as.numeric(as.character(ATACSeq_Gneg1$Abundance))
head(ATACSeq_Gneg1)

#############################
###         Gneg2         ###
#############################

ATACSeq_Gneg2 <- read.delim("ATACSeq_Gneg2.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gneg2 <- ATACSeq_Gneg2[,c(-4,-5,-6)]
colnames(ATACSeq_Gneg2) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gneg2$Abundance <- as.numeric(as.character(ATACSeq_Gneg2$Abundance))
head(ATACSeq_Gneg2); tail(ATACSeq_Gneg2)

#############################
###         Gneg3         ###
#############################

ATACSeq_Gneg3 <- read.delim("ATACSeq_Gneg3.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gneg3 <- ATACSeq_Gneg3[,c(-4,-5,-6)]
colnames(ATACSeq_Gneg3) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gneg3$Abundance <- as.numeric(as.character(ATACSeq_Gneg3$Abundance))
head(ATACSeq_Gneg3); tail(ATACSeq_Gneg3)

#############################
###         Gneg4         ###
#############################

ATACSeq_Gneg4 <- read.delim("ATACSeq_Gneg4.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gneg4 <- ATACSeq_Gneg4[,c(-4,-5,-6)]
colnames(ATACSeq_Gneg4) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gneg4$Abundance <- as.numeric(as.character(ATACSeq_Gneg4$Abundance))
head(ATACSeq_Gneg4); tail(ATACSeq_Gneg4)

#############################
###         Gpos25        ###
#############################

ATACSeq_Gpos25 <- read.delim("ATACSeq_Gpos25.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gpos25 <- ATACSeq_Gpos25[,c(-4,-5)]
colnames(ATACSeq_Gpos25) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gpos25$Abundance <- as.numeric(as.character(ATACSeq_Gpos25$Abundance))
head(ATACSeq_Gpos25)

#############################
###         Gpos50        ###
#############################

ATACSeq_Gpos50 <- read.delim("ATACSeq_Gpos50.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gpos50 <- ATACSeq_Gpos50[,c(-4,-5)]
colnames(ATACSeq_Gpos50) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gpos50$Abundance <- as.numeric(as.character(ATACSeq_Gpos50$Abundance))
head(ATACSeq_Gpos50); tail(ATACSeq_Gpos50)


#############################
###         Gpos75        ###
#############################

ATACSeq_Gpos75 <- read.delim("ATACSeq_Gpos75.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gpos75 <- ATACSeq_Gpos75[,c(-4,-5)]
colnames(ATACSeq_Gpos75) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gpos75$Abundance <- as.numeric(as.character(ATACSeq_Gpos75$Abundance))
head(ATACSeq_Gpos75)

#############################
###         Gpos100       ###
#############################

ATACSeq_Gpos100 <- read.delim("ATACSeq_Gpos100.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gpos100 <- ATACSeq_Gpos100[,c(-4,-5)]
colnames(ATACSeq_Gpos100) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gpos100$Abundance <- as.numeric(as.character(ATACSeq_Gpos100$Abundance))
head(ATACSeq_Gpos100)

setwd("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/2_Abundance_H1s_GenomicStructures/2_output_bedtoolsMap_H1s_GenomicStructures")

# Files obtained from the abundance of H1.2 and H1X in all genomic structures

#############################
###         Gneg1         ###
#############################

H1.2_Gneg1 <- read.delim("H1.2_Gneg1.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_Gneg1 <- H1.2_Gneg1[,c(-4,-5,-6)]
colnames(H1.2_Gneg1) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_Gneg1$Abundance <- as.numeric(as.character(H1.2_Gneg1$Abundance))
head(H1.2_Gneg1)


H1X_Gneg1 <- read.delim("H1X_Gneg1.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_Gneg1 <- H1X_Gneg1[,c(-4,-5,-6)]
colnames(H1X_Gneg1) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_Gneg1$Abundance <- as.numeric(as.character(H1X_Gneg1$Abundance))
head(H1X_Gneg1)

#############################
###         Gneg2         ###
#############################

H1.2_Gneg2 <- read.delim("H1.2_Gneg2.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_Gneg2 <- H1.2_Gneg2[,c(-4,-5,-6)]
colnames(H1.2_Gneg2) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_Gneg2$Abundance <- as.numeric(as.character(H1.2_Gneg2$Abundance))
head(H1.2_Gneg2)

H1X_Gneg2 <- read.delim("H1X_Gneg2.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_Gneg2 <- H1X_Gneg2[,c(-4,-5,-6)]
colnames(H1X_Gneg2) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_Gneg2$Abundance <- as.numeric(as.character(H1X_Gneg2$Abundance))
head(H1X_Gneg2)


#############################
###         Gneg3         ###
#############################

H1.2_Gneg3 <- read.delim("H1.2_Gneg3.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_Gneg3 <- H1.2_Gneg3[,c(-4,-5,-6)]
colnames(H1.2_Gneg3) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_Gneg3$Abundance <- as.numeric(as.character(H1.2_Gneg3$Abundance))
head(H1.2_Gneg3)

H1X_Gneg3 <- read.delim("H1X_Gneg3.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_Gneg3 <- H1X_Gneg3[,c(-4,-5,-6)]
colnames(H1X_Gneg3) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_Gneg3$Abundance <- as.numeric(as.character(H1X_Gneg3$Abundance))
head(H1X_Gneg3)

#############################
###         Gneg4         ###
#############################

H1.2_Gneg4 <- read.delim("H1.2_Gneg4.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_Gneg4 <- H1.2_Gneg4[,c(-4,-5,-6)]
colnames(H1.2_Gneg4) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_Gneg4$Abundance <- as.numeric(as.character(H1.2_Gneg4$Abundance))
head(H1.2_Gneg4)

H1X_Gneg4 <- read.delim("H1X_Gneg4.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_Gneg4 <- H1X_Gneg4[,c(-4,-5,-6)]
colnames(H1X_Gneg4) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_Gneg4$Abundance <- as.numeric(as.character(H1X_Gneg4$Abundance))
head(H1X_Gneg4)


#############################
###         Gpos25        ###
#############################

H1.2_Gpos25 <- read.delim("H1.2_Gpos25.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_Gpos25 <- H1.2_Gpos25[,c(-4,-5)]
colnames(H1.2_Gpos25) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_Gpos25$Abundance <- as.numeric(as.character(H1.2_Gpos25$Abundance))
head(H1.2_Gpos25)

H1X_Gpos25 <- read.delim("H1X_Gpos25.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_Gpos25 <- H1X_Gpos25[,c(-4,-5)]
colnames(H1X_Gpos25) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_Gpos25$Abundance <- as.numeric(as.character(H1X_Gpos25$Abundance))
head(H1X_Gpos25)


#############################
###         Gpos50        ###
#############################

H1.2_Gpos50 <- read.delim("H1.2_Gpos50.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_Gpos50 <- H1.2_Gpos50[,c(-4,-5)]
colnames(H1.2_Gpos50) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_Gpos50$Abundance <- as.numeric(as.character(H1.2_Gpos50$Abundance))
head(H1.2_Gpos50)

H1X_Gpos50 <- read.delim("H1X_Gpos50.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_Gpos50 <- H1X_Gpos50[,c(-4,-5)]
colnames(H1X_Gpos50) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_Gpos50$Abundance <- as.numeric(as.character(H1X_Gpos50$Abundance))
head(H1X_Gpos50)

#############################
###         Gpos75        ###
#############################

H1.2_Gpos75 <- read.delim("H1.2_Gpos75.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_Gpos75 <- H1.2_Gpos75[,c(-4,-5)]
colnames(H1.2_Gpos75) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_Gpos75$Abundance <- as.numeric(as.character(H1.2_Gpos75$Abundance))
head(H1.2_Gpos75)

H1X_Gpos75 <- read.delim("H1X_Gpos75.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_Gpos75 <- H1X_Gpos75[,c(-4,-5)]
colnames(H1X_Gpos75) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_Gpos75$Abundance <- as.numeric(as.character(H1X_Gpos75$Abundance))
head(H1X_Gpos75)

#############################
###         Gpos100       ###
#############################

H1.2_Gpos100 <- read.delim("H1.2_Gpos100.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_Gpos100 <- H1.2_Gpos100[,c(-4,-5)]
colnames(H1.2_Gpos100) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_Gpos100$Abundance <- as.numeric(as.character(H1.2_Gpos100$Abundance))
head(H1.2_Gpos100)

H1X_Gpos100 <- read.delim("H1X_Gpos100.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_Gpos100 <- H1X_Gpos100[,c(-4,-5)]
colnames(H1X_Gpos100) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_Gpos100$Abundance <- as.numeric(as.character(H1X_Gpos100$Abundance))
head(H1X_Gpos100)


#############################
#     Correlation Tests     #    
#############################

# Grouping all G Bands in the same column

Cytoband_Group <- c(rep("Gneg1", 104),rep("Gneg2", 104),rep("Gneg3", 103),rep("Gneg4", 103),rep("Gpos25", 87), rep("Gpos50", 121), rep("Gpos75", 89), rep("Gpos100", 81))
H1.2 <- c(H1.2_Gneg1$Abundance,H1.2_Gneg2$Abundance,H1.2_Gneg3$Abundance,H1.2_Gneg4$Abundance,H1.2_Gpos25$Abundance,H1.2_Gpos50$Abundance,H1.2_Gpos75$Abundance,H1.2_Gpos100$Abundance)
H1X <- c(H1X_Gneg1$Abundance,H1X_Gneg2$Abundance,H1X_Gneg3$Abundance,H1X_Gneg4$Abundance,H1X_Gpos25$Abundance,H1X_Gpos50$Abundance,H1X_Gpos75$Abundance,H1X_Gpos100$Abundance)
ATACseq <- c(ATACSeq_Gneg1$Abundance,ATACSeq_Gneg2$Abundance,ATACSeq_Gneg3$Abundance,ATACSeq_Gneg4$Abundance,ATACSeq_Gpos25$Abundance,ATACSeq_Gpos50$Abundance,ATACSeq_Gpos75$Abundance,ATACSeq_Gpos100$Abundance)

histone1 <- data.frame(cbind(H1.2,H1X,ATACseq,Cytoband_Group))
histone1$H1.2 <- as.numeric(paste0(histone1$H1.2))
histone1$H1X <-  as.numeric(paste0(histone1$H1X))
histone1$ATACseq <-  as.numeric(paste0(histone1$ATACseq))
histone1 <- histone1[which(histone1$ATACseq <= 0.2),]
histone1 <- histone1[which(histone1$H1.2 >= -0.01),]
head(histone1)

# Correlation test with Pearson's method

a <- cor.test(histone1$H1.2, histone1$ATACseq, method= "pearson", na.action = "na.exclude")
b <- cor.test(histone1$H1X, histone1$ATACseq, method= "pearson", na.action = "na.exclude")
  
Correlation_Table1 <- data.frame(x= c(a$p.value,b$p.value),y= c(a$estimate,b$estimate))
Correlation_Table1 <- as.data.frame(t(Correlation_Table1))
row.names(Correlation_Table1) <- c("p.value", "corr.coef")
colnames(Correlation_Table1) <- c("H1.2-ATACseq","H1X-ATACseq")
Correlation_Table1


#####################
#   Scatter Plots   #      
#####################
library(ggplot2)
library(gridExtra)

GBand_H1.2 <- ggplot(histone1, aes(x=H1.2, y=ATACseq)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=Cytoband_Group),  show.legend = T, size=1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "H1.2 vs ATAC-seq", subtitle = "P-value = 4.101944e-44 \nCorrelation coefficient = -4.697294e-01") +
  labs(x = "H1.2",y = "ATAC-seq")

GBand_H1X <- ggplot(histone1, aes(x=H1X, y=ATACseq)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=Cytoband_Group),  show.legend = T, size=1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "H1X vs ATAC-seq", subtitle = "P-value = 1.803749e-131 \nCorrelation coefficient = 7.310765e-01") +
  labs(x = "H1X", y = "ATAC-seq")

grid.arrange(GBand_H1.2,GBand_H1X,nrow=1)



setwd("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/1_ATACseq/2_output_abundance")

# READING BED FILES OF A/B COMPARTMENTS

#############################
###     A compartment     ###
#############################

ATACSeq_A <- read.delim("ATACSeq_A.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_A <- ATACSeq_A[,-4]
colnames(ATACSeq_A) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_A$Abundance <- as.numeric(as.character(ATACSeq_A$Abundance))
head(ATACSeq_A)


#############################
###     B compartment     ###
#############################

ATACSeq_B <- read.delim("ATACSeq_B.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_B <- ATACSeq_B[,-4]
colnames(ATACSeq_B) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_B$Abundance <- as.numeric(as.character(ATACSeq_B$Abundance))
head(ATACSeq_B)


setwd("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/2_Abundance_H1s_GenomicStructures/2_output_bedtoolsMap_H1s_GenomicStructures")

#############################
###     A compartment     ###
#############################

H1.2_A <- read.delim("H1.2_A.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_A <- H1.2_A[,-4]
colnames(H1.2_A) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_A$Abundance <- as.numeric(as.character(H1.2_A$Abundance))
head(H1.2_A)

H1X_A <- read.delim("H1X_A.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_A <- H1X_A[,-4]
colnames(H1X_A) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_A$Abundance <- as.numeric(as.character(H1X_A$Abundance))
head(H1X_A)

#############################
###     B compartment     ###
#############################

H1.2_B <- read.delim("H1.2_B.bed", header=FALSE,stringsAsFactors=FALSE)
H1.2_B <- H1.2_B[,-4]
colnames(H1.2_B) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_B$Abundance <- as.numeric(as.character(H1.2_B$Abundance))
head(H1.2_B)

H1X_B <- read.delim("H1X_B.bed", header=FALSE,stringsAsFactors=FALSE)
H1X_B <- H1X_B[,-4]
colnames(H1X_B) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_B$Abundance <- as.numeric(as.character(H1X_B$Abundance))
head(H1X_B)


#############################
#     Correlation Tests     #    
#############################

# grouping A/B compartments in the same column

Compartment <- c(rep("A", 1098),rep("B", 1098))
H1.2 <- c(H1.2_A$Abundance,H1.2_B$Abundance)
H1X <- c(H1X_A$Abundance,H1X_B$Abundance)
ATACseq <- c(ATACSeq_A$Abundance,ATACSeq_B$Abundance)

histone2 <- data.frame(cbind(H1.2,H1X,ATACseq,Compartment))
histone2$H1.2 <- as.numeric(paste0(histone2$H1.2))
histone2$H1X <-  as.numeric(paste0(histone2$H1X))
histone2$ATACseq <-  as.numeric(paste0(histone2$ATACseq))
histone2 <- histone2[which(histone2$H1.2 >= -0.01),]
histone2 <- histone2[which(histone2$H1.2 <= 0.02),]
histone2 <- histone2[which(histone2$ATACseq <= 0.125),]
head(histone2)

# correlation test with Pearson's method

c <- cor.test(histone2$H1.2, histone2$ATACseq, method= "pearson", na.action = "na.exclude")
d <- cor.test(histone2$H1X, histone2$ATACseq, method= "pearson", na.action = "na.exclude")

Correlation_Table2 <- data.frame(x= c(c$p.value,d$p.value),y= c(c$estimate,d$estimate))
Correlation_Table2 <- as.data.frame(t(Correlation_Table2))
row.names(Correlation_Table2) <- c("p.value", "corr.coef")
colnames(Correlation_Table2) <- c("H1.2-ATACseq","H1X-ATACseq")
Correlation_Table2

#####################
#   Scatter Plots   #      
#####################

ABcomp_H1.2 <- ggplot(histone2, aes(x=H1.2, y=ATACseq)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=Compartment),  show.legend = T, size=1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "H1.2 vs ATAC-seq", subtitle = "P-value = 1.352711e-50 \nCorrelation coefficient = -3.123793e-01") +
  labs(x = "H1.2",y = "ATAC-seq")

ABcomp_H1X <- ggplot(histone2, aes(x=H1X, y=ATACseq)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=Compartment),  show.legend = T, size=1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "H1X vs ATAC-seq", subtitle = "P-value = 5.078062e-140 \nCorrelation coefficient = 5.025902e-01") +
  labs(x = "H1X",y = "ATAC-seq")

grid.arrange(ABcomp_H1.2,ABcomp_H1X,nrow=1)




setwd("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/2_Abundance_H1s_GenomicStructures/2_output_bedtoolsMap_H1s_GenomicStructures")

# READING BED FILES OF TADs

#############################
###         H1.2          ###
#############################

H1.2_H1.2 <- read.delim("H1.2_H1.2.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1.2_H1.2) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_H1.2$Abundance <- as.numeric(as.character(H1.2_H1.2$Abundance))
head(H1.2_H1.2)

H1X_H1.2 <- read.delim("H1X_H1.2.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1X_H1.2) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_H1.2$Abundance <- as.numeric(as.character(H1X_H1.2$Abundance))
head(H1X_H1.2)

#############################
###          H1X          ###
#############################

H1.2_H1X <- read.delim("H1.2_H1X.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1.2_H1X) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_H1X$Abundance <- as.numeric(as.character(H1.2_H1X$Abundance))
head(H1.2_H1X)

H1X_H1X <- read.delim("H1X_H1X.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1X_H1X) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_H1X$Abundance <- as.numeric(as.character(H1X_H1X$Abundance))
head(H1X_H1X)


#############################
###      LIKE_H1.2        ###
#############################

H1.2_LIKE_H1.2 <- read.delim("H1.2_LIKE_H1.2.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1.2_LIKE_H1.2) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_LIKE_H1.2$Abundance <- as.numeric(as.character(H1.2_LIKE_H1.2$Abundance))
head(H1.2_LIKE_H1.2)

H1X_LIKE_H1.2 <- read.delim("H1X_LIKE_H1.2.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1X_LIKE_H1.2) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_LIKE_H1.2$Abundance <- as.numeric(as.character(H1X_LIKE_H1.2$Abundance))
head(H1X_LIKE_H1.2)

#############################
###       LIKE_H1X        ###
#############################

H1.2_LIKE_H1X <- read.delim("H1.2_LIKE_H1X.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1.2_LIKE_H1X) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_LIKE_H1X$Abundance <- as.numeric(as.character(H1.2_LIKE_H1X$Abundance))
head(H1.2_LIKE_H1X)

H1X_LIKE_H1X <- read.delim("H1X_LIKE_H1X.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1X_LIKE_H1X) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_LIKE_H1X$Abundance <- as.numeric(as.character(H1X_LIKE_H1X$Abundance))
head(H1X_LIKE_H1X)

setwd("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/1_ATACseq/2_output_abundance")

#############################
###         H1.2          ###
#############################

ATACSeq_H1.2 <- read.delim("ATACSeq_H1.2.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(ATACSeq_H1.2) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_H1.2$Abundance <- as.numeric(as.character(ATACSeq_H1.2$Abundance))
head(ATACSeq_H1.2)

#############################
###          H1X          ###
#############################

ATACSeq_H1X <- read.delim("ATACSeq_H1X.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(ATACSeq_H1X) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_H1X$Abundance <- as.numeric(as.character(ATACSeq_H1X$Abundance))
head(ATACSeq_H1X)


#############################
###      LIKE_H1.2        ###
#############################

ATACSeq_LIKE_H1.2 <- read.delim("ATACSeq_LIKE_H1.2.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(ATACSeq_LIKE_H1.2) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_LIKE_H1.2$Abundance <- as.numeric(as.character(ATACSeq_LIKE_H1.2$Abundance))
head(ATACSeq_LIKE_H1.2)

#############################
###       LIKE_H1X        ###
#############################

ATACSeq_LIKE_H1X <- read.delim("ATACSeq_LIKE_H1X.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(ATACSeq_LIKE_H1X) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_LIKE_H1X$Abundance <- as.numeric(as.character(ATACSeq_LIKE_H1X$Abundance))
head(ATACSeq_LIKE_H1X)


#############################
#     Correlation Tests     #    
#############################

# grouping all TADs in the same column

TADs <- c(rep("H1.2", 802),rep("H1X", 803),rep("LIKE_H1.2", 802),rep("LIKE_H1X", 802))
H1.2 <- c(H1.2_H1.2$Abundance,H1.2_H1X$Abundance,H1.2_LIKE_H1.2$Abundance,H1.2_LIKE_H1X$Abundance)
H1X <- c(H1X_H1.2$Abundance,H1X_H1X$Abundance,H1X_LIKE_H1.2$Abundance,H1X_LIKE_H1X$Abundance)
ATACseq <- c(ATACSeq_H1.2$Abundance,ATACSeq_H1X$Abundance,ATACSeq_LIKE_H1.2$Abundance,ATACSeq_LIKE_H1X$Abundance)

histone3 <- data.frame(cbind(H1.2,H1X,ATACseq,TADs))
histone3$H1.2 <- as.numeric(paste0(histone3$H1.2))
histone3$H1X <-  as.numeric(paste0(histone3$H1X))
histone3$ATACseq <-  as.numeric(paste0(histone3$ATACseq))
histone3 <- histone3[which(histone3$H1.2 >= -0.0125),]
histone3 <- histone3[which(histone3$ATACseq <= 0.125),]
histone3 <- histone3[which(histone3$H1X <= 0.02),]
head(histone3)

# correlation test with Pearson's method

e <- cor.test(histone3$H1.2, histone3$ATACseq, method= "pearson")
f <- cor.test(histone3$H1X, histone3$ATACseq, method= "pearson")

Correlation_Table3 <- data.frame(x= c(e$p.value,f$p.value),y= c(e$estimate,f$estimate))
Correlation_Table3 <- as.data.frame(t(Correlation_Table3))
row.names(Correlation_Table3) <- c("p.value", "corr.coef")
colnames(Correlation_Table3) <- c("H1.2-ATACseq","H1X-ATACseq")
Correlation_Table3


#####################
#   Scatter Plots   #      
#####################

TADs_H1.2 <- ggplot(histone3, aes(x=H1.2, y=ATACseq)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=TADs),  show.legend = T, size=1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "H1.2 vs ATACseq", subtitle = "P-value = 3.490877e-164 \nCorrelation coefficient = -4.561442e-01") +
  labs(x = "H1.2",y = "ATAC-seq")

TADs_H1X <- ggplot(histone3, aes(x=H1X, y=ATACseq)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=TADs),  show.legend = T, size=1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "H1X vs ATACseq", subtitle = "P-value = 0.0000000 \nCorrelation coefficient = 0.6828587") +
  labs(x = "H1X",y = "ATAC-seq")

grid.arrange(TADs_H1.2,TADs_H1X,nrow=1)
