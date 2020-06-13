setwd("~/Documents/Bioinformatica/4🍟/TFG/Retrovirus/4_Abundance_H1s_Virus/3_output_bedtoolsMap")

H1.2_HIV <- read.delim("H1.2_HIV.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1.2_HIV) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_HIV$Abundance <- as.numeric(as.character(H1.2_HIV$Abundance))
head(H1.2_HIV)

H1.2_HTLV <- read.delim("H1.2_HTLV.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1.2_HTLV) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1.2_HTLV$Abundance <- as.numeric(as.character(H1.2_HTLV$Abundance))
head(H1.2_HTLV)

H1X_HIV <- read.delim("H1X_HIV.bed", header=FALSE,stringsAsFactors=FALSE)
colnames(H1X_HIV) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_HIV$Abundance <- as.numeric(as.character(H1X_HIV$Abundance))
head(H1X_HIV)

H1X_HTLV <- read.delim("H1X_HTLV.bed", header=FALSE, stringsAsFactors=FALSE)
colnames(H1X_HTLV) <- c("Chromosome","StartCoord","EndCoord","Abundance")
H1X_HTLV$Abundance <- as.numeric(as.character(H1X_HTLV$Abundance))
head(H1X_HTLV)


############ PLOT ABUNDANCE OF H1.2 AND H1X IN HIV ############

Colors <- c("lightblue2","sandybrown")

boxplot(H1.2_HIV$Abundance,H1X_HIV$Abundance,
        col = Colors, outline = F, main = "Abundance of H1.2 and H1X in HIV", 
        ylab="Abundance",  at=c(1,3), xaxt="n")

abline(v=2, col = "grey", lwd = 1)

title(sub="H1.2", adj=0.2, line=0.8)
title(sub="H1X", adj=0.85, line=0.8)

############ PLOT ABUNDANCE OF H1.2 AND H1X IN HTLV ############

Colors <- c("lightblue2","sandybrown")

boxplot(H1.2_HTLV$Abundance,H1X_HTLV$Abundance,
        col = Colors, outline = F, main = "Abundance of H1.2 and H1X in HTLV", 
        ylab="Abundance",  at=c(1,3), xaxt="n")

abline(v=2, col = "grey", lwd = 1)

title(sub="H1.2", adj=0.2, line=0.8)
title(sub="H1X", adj=0.85, line=0.8)
