setwd("~/Documents/Bioinformatica/4🍟/TFG/Retrovirus/2_Abundance_H1s_GenomicStructures/2_output_bedtoolsMap_H1s_GenomicStructures")
library(ggplot2)

# READ BED FILES OF G BANDS

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

# grouping all G Bands in the same column

Cytoband_Group <- c(rep("Gneg1", 104),rep("Gneg2", 104),rep("Gneg3", 103),rep("Gneg4", 103),rep("Gpos25", 87), rep("Gpos50", 121), rep("Gpos75", 89), rep("Gpos100", 81))
H1.2 <- c(H1.2_Gneg1$Abundance,H1.2_Gneg2$Abundance,H1.2_Gneg3$Abundance,H1.2_Gneg4$Abundance,H1.2_Gpos25$Abundance,H1.2_Gpos50$Abundance,H1.2_Gpos75$Abundance,H1.2_Gpos100$Abundance)
H1X <- c(H1X_Gneg1$Abundance,H1X_Gneg2$Abundance,H1X_Gneg3$Abundance,H1X_Gneg4$Abundance,H1X_Gpos25$Abundance,H1X_Gpos50$Abundance,H1X_Gpos75$Abundance,H1X_Gpos100$Abundance)

histone1 <- data.frame(cbind(H1.2,H1X,Cytoband_Group))
histone1$H1.2 <- as.numeric(paste0(histone1$H1.2))
histone1$H1X <-  as.numeric(paste0(histone1$H1X))
histone1 <- histone1[which(histone1$H1.2 >= -0.01),]
head(histone1)

# correlation test with Pearson's method

a <- cor.test(histone1$H1.2, histone1$H1X, method= "pearson")

Correlation_Table <- data.frame(x= c(a$p.value),y= c(a$estimate))
Correlation_Table <- as.data.frame(t(Correlation_Table))
row.names(Correlation_Table) <- c("p.value", "corr.coef")
Correlation_Table

#####################
#   Scatter Plots   #      
#####################

ggplot(histone1, aes(x=H1.2, y=H1X)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=Cytoband_Group),  show.legend = T, size=1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "H1.2 vs H1X", subtitle = "P-value = 1.141982e-41 \nCorrelation coefficient = -4.554800e-01") +
  labs(x="H1.2",y="H1X")

# READ BED FILES OF A/B COMPARTMENTS

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

histone2 <- data.frame(cbind(H1.2,H1X,Compartment))
histone2$H1.2 <- as.numeric(paste0(histone2$H1.2))
histone2$H1X <-  as.numeric(paste0(histone2$H1X))
histone2 <- histone2[which(histone2$H1.2 >= -0.005),]
histone2 <- histone2[which(histone2$H1X <= 0.01),]
histone2 <- histone2[which(histone2$H1X >= -0.01),]
head(histone2)

# correlation test with Pearson's method

a <- cor.test(histone2$H1.2, histone2$H1X, method= "pearson")

Correlation_Table <- data.frame(x= c(a$p.value),y= c(a$estimate))
Correlation_Table <- as.data.frame(t(Correlation_Table))
row.names(Correlation_Table) <- c("p.value", "corr.coef")
Correlation_Table

#####################
#   Scatter Plots   #      
#####################

ggplot(histone2, aes(x=H1.2, y=H1X)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=Compartment),  show.legend = T, size=1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "H1.2 vs H1X", subtitle = "P-value = 6.633118e-07 \nCorrelation coefficient = -1.070784e-01") +
  labs(x="H1.2",y="H1X")

# READ BED FILES OF TADs

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

#############################
#     Correlation Tests     #    
#############################

# grouping all TADs in the same column

TADs <- c(rep("H1.2", 802),rep("H1X", 803),rep("LIKE_H1.2", 802),rep("LIKE_H1X", 802))
H1.2 <- c(H1.2_H1.2$Abundance,H1.2_H1X$Abundance,H1.2_LIKE_H1.2$Abundance,H1.2_LIKE_H1X$Abundance)
H1X <- c(H1X_H1.2$Abundance,H1X_H1X$Abundance,H1X_LIKE_H1.2$Abundance,H1X_LIKE_H1X$Abundance)

histone3 <- data.frame(cbind(H1.2,H1X,TADs))
histone3$H1.2 <- as.numeric(paste0(histone3$H1.2))
histone3$H1X <-  as.numeric(paste0(histone3$H1X))
histone3 <- histone3[which(histone3$H1.2 >= -0.003),]
histone3 <- histone3[which(histone3$H1.2 <= 0.0051),]
histone3 <- histone3[which(histone3$H1X <= 0.008),]
histone3 <- histone3[which(histone3$H1X > -0.01),]
head(histone3)

# correlation test with Pearson's method

a <- cor.test(histone3$H1.2, histone3$H1X, method= "pearson")

Correlation_Table <- data.frame(x= c(a$p.value),y= c(a$estimate))
Correlation_Table <- as.data.frame(t(Correlation_Table))
row.names(Correlation_Table) <- c("p.value", "corr.coef")
Correlation_Table

#####################
#   Scatter Plots   #      
#####################

ggplot(histone3, aes(x=H1.2, y=H1X)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=TADs),  show.legend = T, size=1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "H1.2 vs H1X", subtitle = "P-value = 6.589722e-67 \nCorrelation coefficient = -3.107523e-01") +
  labs(x="H1.2",y="H1X")
