setwd("~/Desktop")
library(readr)

# creating bed files of HIV and HTLV

#############################
###          HIV          ###
#############################

chr_HIV <- read.table("chr_HIV.txt", header = T)
colnames(chr_HIV) <- "Chromosome"
head(chr_HIV)

insertPos_HIV <- read.table("insertPos_HIV.txt", header = T)
colnames(insertPos_HIV) <- "Insert_Position"
head(insertPos_HIV)

HIV <- data.frame(cbind(chr_HIV,insertPos_HIV))
HIV$Insert_Position <- as.numeric(as.character(HIV$Insert_Position))
head(HIV)


StartCoord <- c()
EndCoord <- c()
Chromosome <- c()
start <- 0
end <- 0

for (i in seq(nrow(HIV))){
  start <- HIV[i,2] - 150 # we subtract 150 bp to the insertion position to obtain the start coordinate
  StartCoord <- c(StartCoord,start)
  start <- 0
  end <- HIV[i,2] + 150 # we add 150 bp to the insertion position to obtain the end coordinate
  EndCoord <- c(EndCoord,end)
  end <- 0
  Chromosome <- c(Chromosome, HIV[i,1])
}

HIV_bed <- data.frame(cbind(Chromosome,StartCoord,EndCoord))
HIV_bed <- HIV_bed[order(HIV_bed$Chromosome,HIV_bed$StartCoord),]
head(HIV_bed)

write.table(HIV_bed, "/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/Abundance/HIV_bed.bed", quote = F, sep = "\t", row.names = F, col.names = F)

#############################
###         HTLV          ###
#############################

chr_HTLV <- read.table("chr_HTLV.txt", header = T)
colnames(chr_HTLV) <- "Chromosome"
head(chr_HTLV)

insertPos_HTLV <- read.table("insertPos_HTLV.txt", header = T)
colnames(insertPos_HTLV) <- "Insert_Position"
head(insertPos_HTLV)

HTLV <- data.frame(cbind(chr_HTLV,insertPos_HTLV))
HTLV$Insert_Position <- as.numeric(as.character(HTLV$Insert_Position))
head(HTLV)


StartCoord <- c()
EndCoord <- c()
Chromosome <- as.factor(c())
start <- 0
end <- 0

for (i in seq(nrow(HTLV))){
  start <- HTLV[i,2] - 150 # we subtract 150 bp to the insertion position to obtain the start coordinate
  StartCoord <- c(StartCoord,start)
  start <- 0
  end <- HTLV[i,2] + 150 # we add 150 bp to the insertion position to obtain the end coordinate
  EndCoord <- c(EndCoord,end)
  end <- 0
  Chromosome <- as.factor(c(Chromosome, HTLV[i,1]))
}

HTLV_bed <- data.frame(cbind(Chromosome,StartCoord,EndCoord))
HTLV_bed <- HTLV_bed[order(HTLV_bed$Chromosome,HTLV_bed$StartCoord),]
head(HTLV_bed)

write.table(HTLV_bed, "/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/Abundance/HTLV_bed.bed", quote = F, sep = "\t", row.names = F, col.names = F)
