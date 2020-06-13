setwd("~/Documents/Bioinformatica/4🍟/TFG/Retrovirus/2_Abundance_H1s_GenomicStructures/2_output_bedtoolsMap_H1s_GenomicStructures")

# READING BED FILES OF G BANDS

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

# CREATING BOXPLOT OF ABUNDANCE BETWEEN H1.2 and H1X in G Bands

GBands_Colors <- c("lightblue2","sandybrown","lightblue2","sandybrown","lightblue2","sandybrown","lightblue2","sandybrown","lightblue2","sandybrown","lightblue2","sandybrown","lightblue2","sandybrown","lightblue2","sandybrown")

boxplot(H1.2_Gneg1$Abundance,H1X_Gneg1$Abundance,
        H1.2_Gneg2$Abundance,H1X_Gneg2$Abundance,
        H1.2_Gneg3$Abundance,H1X_Gneg3$Abundance,
        H1.2_Gneg4$Abundance,H1X_Gneg4$Abundance,
        H1.2_Gpos25$Abundance,H1X_Gpos25$Abundance,
        H1.2_Gpos50$Abundance,H1X_Gpos50$Abundance,
        H1.2_Gpos75$Abundance,H1X_Gpos75$Abundance,
        H1.2_Gpos100$Abundance,H1X_Gpos100$Abundance,
        col = GBands_Colors, outline = F, main = "Abundance of H1.2 and H1X in GBands", 
        ylab="Abundance",  at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23), xaxt="n")

abline(v=3, col = "grey", lwd = 1)
abline(v=6, col = "grey", lwd = 1)
abline(v=9, col = "grey", lwd = 1)
abline(v=12, col = "grey", lwd = 1)
abline(v=15, col = "grey", lwd = 1)
abline(v=18, col = "grey", lwd = 1)
abline(v=21, col = "grey", lwd = 1)

title(sub="Gneg1", adj=0.045, line=0.8)
title(sub="Gneg2", adj=0.18, line=0.8)
title(sub="Gneg3", adj=0.31, line=0.8)
title(sub="Gneg4", adj=0.44, line=0.8)
title(sub="Gpos25", adj=0.57, line=0.8)
title(sub="Gpos50", adj=0.7, line=0.8)
title(sub="Gpos75", adj=0.84, line=0.8)
title(sub="Gpos100", adj=0.98, line=0.8)

legend("topright", legend = c("H1.2","H1X"), fill=c("lightblue2","sandybrown"),cex=0.9, inset=0.01)

# READING BED FILES OF A/B COMPARTMENTS

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

# CREATING BOXPLOT OF ABUNDANCE BETWEEN H1.2 and H1X in A/B compartments

AB_Colors <- c("lightblue2","sandybrown","lightblue2","sandybrown")

boxplot(H1.2_A$Abundance,H1X_A$Abundance,
        H1.2_B$Abundance,H1X_B$Abundance,
        col = AB_Colors, outline = F, main = "Abundance of H1.2 and H1X in A/B compartments", 
        ylab="Abundance",  at=c(1,2,4,5), xaxt="n")

abline(v=3, col = "grey", lwd = 1)

title(sub="A compartment", adj=0.2, line=0.8)
title(sub="B compartment", adj=0.85, line=0.8)

legend("topleft", legend = c("H1.2","H1X"), fill=c("lightblue2","sandybrown"),cex=1, inset=0.01)

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

# CREATING BOXPLOT OF ABUNDANCE BETWEEN H1.2 and H1X in TADs

TADs_Colors <- c("lightblue2","sandybrown","lightblue2","sandybrown","lightblue2","sandybrown","lightblue2","sandybrown")

boxplot(H1.2_H1X$Abundance,H1X_H1X$Abundance,
        H1.2_LIKE_H1X$Abundance,H1X_LIKE_H1X$Abundance,
        H1.2_LIKE_H1.2$Abundance,H1X_LIKE_H1.2$Abundance,
        H1.2_H1.2$Abundance,H1X_H1.2$Abundance,
        col = TADs_Colors, outline = F, main = "Abundance of H1.2 and H1X in TADs", 
        ylab="Abundance",  at=c(1,2,4,5,7,8,10,11), xaxt="n")

abline(v=3, col = "grey", lwd = 1)
abline(v=6, col = "grey", lwd = 1)
abline(v=9, col = "grey", lwd = 1)

title(sub="H1X", adj=0.1, line=0.8)
title(sub="LIKE_H1X", adj=0.36, line=0.8)
title(sub="LIKE_H1.2", adj=0.65, line=0.8)
title(sub="H1.2", adj=0.91, line=0.8)

legend("topright", legend = c("H1.2","H1X"), fill=c("lightblue2","sandybrown"),cex=0.9, inset=0.01)

