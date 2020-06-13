setwd("~/Documents/Bioinformatica/4🍟/TFG/Retrovirus/1_ATACseq/2_output_abundance")

# READING BED FILES

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
ATACSeq_Gneg2 <- ATACSeq_Gneg2[-104,] # removing rows with NAs
head(ATACSeq_Gneg2); tail(ATACSeq_Gneg2)

#############################
###         Gneg3         ###
#############################

ATACSeq_Gneg3 <- read.delim("ATACSeq_Gneg3.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gneg3 <- ATACSeq_Gneg3[,c(-4,-5,-6)]
colnames(ATACSeq_Gneg3) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gneg3$Abundance <- as.numeric(as.character(ATACSeq_Gneg3$Abundance))
ATACSeq_Gneg3 <- ATACSeq_Gneg3[c(-102,-103),] # removing rows with NAs
head(ATACSeq_Gneg3); tail(ATACSeq_Gneg3)

#############################
###         Gneg4         ###
#############################

ATACSeq_Gneg4 <- read.delim("ATACSeq_Gneg4.bed", header=FALSE,stringsAsFactors=FALSE)
ATACSeq_Gneg4 <- ATACSeq_Gneg4[,c(-4,-5,-6)]
colnames(ATACSeq_Gneg4) <- c("Chromosome","StartCoord","EndCoord","Abundance")
ATACSeq_Gneg4$Abundance <- as.numeric(as.character(ATACSeq_Gneg4$Abundance))
ATACSeq_Gneg4 <- ATACSeq_Gneg4[c(-102,-103),] # removing rows with NAs
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
ATACSeq_Gpos50 <- ATACSeq_Gpos50[c(-119,-120,-121),] # removing rows with NAs
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

#####################
# CREATING BLOXPLOT #
#####################

GBands_Colors <- c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc","Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3")

boxplot(ATACSeq_Gneg1$Abundance,ATACSeq_Gneg2$Abundance,ATACSeq_Gneg3$Abundance,ATACSeq_Gneg4$Abundance,
        ATACSeq_Gpos25$Abundance,ATACSeq_Gpos50$Abundance,ATACSeq_Gpos75$Abundance,ATACSeq_Gpos100$Abundance,
        col = GBands_Colors, outline = F, main = "Abundance of ATAC-seq in GBands", 
        ylab="Abundance",  at=c(1,3,5,7,9,11,13,15), xaxt="n")

abline(v=2, col = "grey", lwd = 1)
abline(v=4, col = "grey", lwd = 1)
abline(v=6, col = "grey", lwd = 1)
abline(v=8, col = "grey", lwd = 1)
abline(v=10, col = "grey", lwd = 1)
abline(v=12, col = "grey", lwd = 1)
abline(v=14, col = "grey", lwd = 1)

title(sub="Gneg1", adj=0.045, line=0.8)
title(sub="Gneg2", adj=0.18, line=0.8)
title(sub="Gneg3", adj=0.31, line=0.8)
title(sub="Gneg4", adj=0.44, line=0.8)
title(sub="Gpos25", adj=0.57, line=0.8)
title(sub="Gpos50", adj=0.7, line=0.8)
title(sub="Gpos75", adj=0.84, line=0.8)
title(sub="Gpos100", adj=0.98, line=0.8)


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

#####################
# CREATING BLOXPLOT #
#####################

AB_Colors <- c("A_compartment"="sandybrown","B_compartment"="lightblue2", "AB_compartment"="peachpuff")

boxplot(ATACSeq_A$Abundance,ATACSeq_B$Abundance,
        col = AB_Colors, outline = F, main = "Abundance of ATAC-seq in A/B compartments", 
        ylab="Abundance",  at=c(1,3), xaxt="n")

abline(v=2, col = "grey", lwd = 1)

title(sub="A compartment", adj=0.2, line=0.8)
title(sub="B compartment", adj=0.85, line=0.8)


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


#####################
# CREATING BLOXPLOT #
#####################

TADs_Colors <- c("H1X_enriched"="sandybrown","LIKE_H1X_enriched"="coral2","H1.2_enriched"="lightblue2","LIKE_H1.2_enriched"="cyan3")

boxplot(ATACSeq_H1X$Abundance,ATACSeq_LIKE_H1X$Abundance,ATACSeq_LIKE_H1.2$Abundance,ATACSeq_H1.2$Abundance,
        col = TADs_Colors, outline = F, main = "Abundance of ATAC-seq in TADs", 
        ylab="Abundance",  at=c(1,3,5,7), xaxt="n")

abline(v=2, col = "grey", lwd = 1)
abline(v=4, col = "grey", lwd = 1)
abline(v=6, col = "grey", lwd = 1)

title(sub="H1X", adj=0.1, line=0.8)
title(sub="LIKE_H1X", adj=0.36, line=0.8)
title(sub="LIKE_H1.2", adj=0.65, line=0.8)
title(sub="H1.2", adj=0.91, line=0.8)

