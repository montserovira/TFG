setwd("~/Desktop")
library(readr)

#############################
###          HIV          ###
#############################

chr_HIV <- read.table("chr_HIV.txt", header = T) # reading chromosomes of HIV
colnames(chr_HIV) <- "Chromosome"
head(chr_HIV)

insertPos_HIV <- read.table("insertPos_HIV.txt", header = T) # reading insert positions of HIV
colnames(insertPos_HIV) <- "Insert_Position"
head(insertPos_HIV)

HIV <- data.frame(cbind(chr_HIV,insertPos_HIV)) # join chromosomes and insert positions in a dataframe
HIV$Insert_Position <- as.numeric(as.character(HIV$Insert_Position))
head(HIV)
save(HIV, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/3_InsertionSizes/4_plots_InsertPos_HIV/HIV.RData")

#############################
###         HTLV          ###
#############################

chr_HTLV <- read.table("chr_HTLV.txt", header = T) # reading chromosomes of HTLV
colnames(chr_HTLV) <- "Chromosome"
head(chr_HTLV)

insertPos_HTLV <- read.table("insertPos_HTLV.txt", header = T) # reading insert positions of HTLV
colnames(insertPos_HTLV) <- "Insert_Position"
head(insertPos_HTLV)

HTLV <- data.frame(cbind(chr_HTLV,insertPos_HTLV)) # join chromosomes and insert positions in a dataframe
HTLV$Insert_Position <- as.numeric(as.character(HTLV$Insert_Position))
head(HTLV)
save(HTLV, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/3_InsertionSizes/5_plots_InsertPos_HTLV/HTLV.RData")

# READING BED FILES OF ALL GENOMIC STRUCTURES

setwd("~/Documents/Bioinformatica/4🍟/TFG/Retrovirus/Genomic Structures/GBands")

#############################
###         Gneg1         ###
#############################

Gneg1 <- read.table("Gneg1.bed",header = F, stringsAsFactors = F)
Gneg1$V7 <- Gneg1$V3 - Gneg1$V2
Gneg1$V8 <- c(rep("Gneg1",104))
Gneg1 <- Gneg1[,-5]
colnames(Gneg1) <- c("Chromosome","StartCoord", "EndCoord", "Band", "Region", "Difference","Type")
head(Gneg1)

#############################
###         Gneg2         ###
#############################

Gneg2 <- read.table("Gneg2.bed",header = F, stringsAsFactors = F)
Gneg2$V7 <- Gneg2$V3 - Gneg2$V2
Gneg2$V8 <- c(rep("Gneg2",104))
Gneg2 <- Gneg2[,-5]
colnames(Gneg2) <- c("Chromosome","StartCoord", "EndCoord", "Band", "Region", "Difference","Type")
head(Gneg2)

#############################
###         Gneg3         ###
#############################

Gneg3 <- read.table("Gneg3.bed",header = F, stringsAsFactors = F)
Gneg3$V7 <- Gneg3$V3 - Gneg3$V2
Gneg3$V8 <- c(rep("Gneg3",103))
Gneg3 <- Gneg3[,-5]
colnames(Gneg3) <- c("Chromosome","StartCoord", "EndCoord", "Band", "Region", "Difference","Type")
head(Gneg3)

#############################
###         Gneg4         ###
#############################

Gneg4 <- read.table("Gneg4.bed",header = F, stringsAsFactors = F)
Gneg4$V7 <- Gneg4$V3 - Gneg4$V2
Gneg4$V8 <- c(rep("Gneg4",103))
Gneg4 <- Gneg4[,-5]
colnames(Gneg4) <- c("Chromosome","StartCoord", "EndCoord", "Band", "Region", "Difference","Type")
head(Gneg4)

#############################
###         Gpos25        ###
#############################

Gpos25 <- read.table("Gpos25.bed",header = F, stringsAsFactors = F)
Gpos25$V6 <- Gpos25$V3 - Gpos25$V2
Gpos25$V7 <- c(rep("Gpos25",87))
Gpos25 <- Gpos25[,-5]
colnames(Gpos25) <- c("Chromosome","StartCoord", "EndCoord", "Band", "Difference","Type")
head(Gpos25)

#############################
###         Gpos50        ###
#############################

Gpos50 <- read.table("Gpos50.bed",header = F, stringsAsFactors = F)
Gpos50$V6 <- Gpos50$V3 - Gpos50$V2
Gpos50$V7 <- c(rep("Gpos50",121))
Gpos50 <- Gpos50[,-5]
colnames(Gpos50) <- c("Chromosome","StartCoord", "EndCoord", "Band", "Difference","Type")
head(Gpos50)

#############################
###         Gpos75        ###
#############################

Gpos75 <- read.table("Gpos75.bed",header = F, stringsAsFactors = F)
Gpos75$V6 <- Gpos75$V3 - Gpos75$V2
Gpos75$V7 <- c(rep("Gpos75",89))
Gpos75 <- Gpos75[,-5]
colnames(Gpos75) <- c("Chromosome","StartCoord", "EndCoord", "Band", "Difference","Type")
head(Gpos75)

#############################
###         Gpos100       ###
#############################

Gpos100 <- read.table("Gpos100.bed",header = F, stringsAsFactors = F)
Gpos100$V6 <- Gpos100$V3 - Gpos100$V2
Gpos100$V7 <- c(rep("Gpos100",81))
Gpos100 <- Gpos100[,-5]
colnames(Gpos100) <- c("Chromosome","StartCoord", "EndCoord", "Band", "Difference","Type")
head(Gpos100)

setwd("~/Documents/Bioinformatica/4🍟/TFG/Retrovirus/Genomic Structures/AB compartments")

#############################
###     Acompartment      ###
#############################

A <- read.table("Acompartment.bed",header = F, stringsAsFactors = F)
A$V5 <- A$V3 - A$V2
colnames(A) <- c("Chromosome","StartCoord", "EndCoord", "Type","Difference")
head(A)

#############################
###     Bcompartment      ###
#############################

B <- read.table("Bcompartment.bed",header = F, stringsAsFactors = F)
B$V5 <- B$V3 - B$V2
colnames(B) <- c("Chromosome","StartCoord", "EndCoord", "Type","Difference")
head(B)

setwd("~/Documents/Bioinformatica/4🍟/TFG/Retrovirus/Genomic Structures/TADs")

#############################
###         H1.2          ###
#############################

H1.2 <- read.table("H1.2_enriched.bed",header = F, stringsAsFactors = F)
H1.2$V4 <- H1.2$V3 - H1.2$V2
H1.2$V5 <- c(rep("H1.2",802))
colnames(H1.2) <- c("Chromosome","StartCoord", "EndCoord", "Difference","Type")
head(H1.2)

#############################
###         H1X          ###
#############################

H1X <- read.table("H1X_enriched.bed",header = F, stringsAsFactors = F)
H1X$V4 <- H1X$V3 - H1X$V2
H1X$V5 <- c(rep("H1X",803))
colnames(H1X) <- c("Chromosome","StartCoord", "EndCoord", "Difference","Type")
head(H1X)

#############################
###       LIKE H1.2       ###
#############################

LIKE_H1.2 <- read.table("LIKE_H1.2_enriched.bed",header = F, stringsAsFactors = F)
LIKE_H1.2$V4 <- LIKE_H1.2$V3 - LIKE_H1.2$V2
LIKE_H1.2$V5 <- c(rep("LIKE_H1.2",802))
colnames(LIKE_H1.2) <- c("Chromosome","StartCoord", "EndCoord", "Difference","Type")
head(LIKE_H1.2)

#############################
###       LIKE H1X       ###
#############################

LIKE_H1X <- read.table("LIKE_H1X_enriched.bed",header = F, stringsAsFactors = F)
LIKE_H1X$V4 <- LIKE_H1X$V3 - LIKE_H1X$V2
LIKE_H1X$V5 <- c(rep("LIKE_H1X",802))
colnames(LIKE_H1X) <- c("Chromosome","StartCoord", "EndCoord", "Difference","Type")
head(LIKE_H1X)

### CALCULATING THE NUMBER OF INSERT POSITION OF EACH RETROVIRUS IN ALL GENOMIC STRUCTURES

#############################
###       HIV_Gneg1       ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(Gneg1))) {
    if (HIV[i,1] == Gneg1[k,1] && HIV[i,2] >= Gneg1[k,2] && HIV[i,2] <= Gneg1[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, Gneg1[k,2])
        EndCoord <- c(EndCoord, Gneg1[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_Gneg1 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_Gneg1 <- HIV_Gneg1[order(HIV_Gneg1$Chromosome,HIV_Gneg1$StartCoord),]
save(HIV_Gneg1, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_Gneg1.RData")
head(HIV_Gneg1)

#############################
###       HIV_Gneg2       ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(Gneg2))) {
    if (HIV[i,1] == Gneg2[k,1] && HIV[i,2] >= Gneg2[k,2] && HIV[i,2] <= Gneg2[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, Gneg2[k,2])
        EndCoord <- c(EndCoord, Gneg2[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_Gneg2 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_Gneg2 <- HIV_Gneg2[order(HIV_Gneg2$Chromosome,HIV_Gneg2$StartCoord),]
head(HIV_Gneg2)
save(HIV_Gneg2, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_Gneg2.RData")



#############################
###       HIV_Gneg3       ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(Gneg3))) {
    if (HIV[i,1] == Gneg3[k,1] && HIV[i,2] >= Gneg3[k,2] && HIV[i,2] <= Gneg3[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, Gneg3[k,2])
        EndCoord <- c(EndCoord, Gneg3[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_Gneg3 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_Gneg3 <- HIV_Gneg3[order(HIV_Gneg3$Chromosome,HIV_Gneg3$StartCoord),]
head(HIV_Gneg3)
save(HIV_Gneg3, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_Gneg3.RData")


#############################
###       HIV_Gneg4       ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(Gneg4))) {
    if (HIV[i,1] == Gneg4[k,1] && HIV[i,2] >= Gneg4[k,2] && HIV[i,2] <= Gneg4[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, Gneg4[k,2])
        EndCoord <- c(EndCoord, Gneg4[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_Gneg4 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_Gneg4 <- HIV_Gneg4[order(HIV_Gneg4$Chromosome,HIV_Gneg4$StartCoord),]
head(HIV_Gneg4)
save(HIV_Gneg4, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_Gneg4.RData")

#############################
###      HIV_Gpos25       ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(Gpos25))) {
    if (HIV[i,1] == Gpos25[k,1] && HIV[i,2] >= Gpos25[k,2] && HIV[i,2] <= Gpos25[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, Gpos25[k,2])
        EndCoord <- c(EndCoord, Gpos25[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_Gpos25 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_Gpos25 <- HIV_Gpos25[order(HIV_Gpos25$Chromosome,HIV_Gpos25$StartCoord),]
head(HIV_Gpos25)
save(HIV_Gpos25, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_Gpos25.RData")


#############################
###      HIV_Gpos50       ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(Gpos50))) {
    if (HIV[i,1] == Gpos50[k,1] && HIV[i,2] >= Gpos50[k,2] && HIV[i,2] <= Gpos50[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, Gpos50[k,2])
        EndCoord <- c(EndCoord, Gpos50[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_Gpos50 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_Gpos50 <- HIV_Gpos50[order(HIV_Gpos50$Chromosome,HIV_Gpos50$StartCoord),]
head(HIV_Gpos50)
save(HIV_Gpos50, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_Gpos50.RData")

#############################
###      HIV_Gpos75       ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(Gpos75))) {
    if (HIV[i,1] == Gpos75[k,1] && HIV[i,2] >= Gpos75[k,2] && HIV[i,2] <= Gpos75[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, Gpos75[k,2])
        EndCoord <- c(EndCoord, Gpos75[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_Gpos75 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_Gpos75 <- HIV_Gpos75[order(HIV_Gpos75$Chromosome,HIV_Gpos75$StartCoord),]
head(HIV_Gpos75)
save(HIV_Gpos75, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_Gpos75.RData")


#############################
###      HIV_Gpos100      ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(Gpos100))) {
    if (HIV[i,1] == Gpos100[k,1] && HIV[i,2] >= Gpos100[k,2] && HIV[i,2] <= Gpos100[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, Gpos100[k,2])
        EndCoord <- c(EndCoord, Gpos100[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_Gpos100 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_Gpos100 <- HIV_Gpos100[order(HIV_Gpos100$Chromosome,HIV_Gpos100$StartCoord),]
head(HIV_Gpos100)
save(HIV_Gpos100, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_Gpos100.RData")


#############################
###         HIV_A         ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(A))) {
    if (HIV[i,1] == A[k,1] && HIV[i,2] >= A[k,2] && HIV[i,2] <= A[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the A compartment
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, A[k,2])
        EndCoord <- c(EndCoord, A[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_A <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_A <- HIV_A[order(HIV_A$Chromosome,HIV_A$StartCoord),]
head(HIV_A)
save(HIV_A, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_A.RData")

#############################
###         HIV_B         ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(B))) {
    if (HIV[i,1] == B[k,1] && HIV[i,2] >= B[k,2] && HIV[i,2] <= B[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the B compartment
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, B[k,2])
        EndCoord <- c(EndCoord, B[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_B <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_B <- HIV_B[order(HIV_B$Chromosome,HIV_B$StartCoord),]
head(HIV_B)
save(HIV_B, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_B.RData")


#############################
###       HIV_H1.2        ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(H1.2))) {
    if (HIV[i,1] == H1.2[k,1] && HIV[i,2] >= H1.2[k,2] && HIV[i,2] <= H1.2[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the TAD
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, H1.2[k,2])
        EndCoord <- c(EndCoord, H1.2[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_H1.2 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_H1.2 <- HIV_H1.2[order(HIV_H1.2$Chromosome,HIV_H1.2$StartCoord),]
head(HIV_H1.2)
save(HIV_H1.2, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_H1.2.RData")

#############################
###       HIV_H1X        ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(H1X))) {
    if (HIV[i,1] == H1X[k,1] && HIV[i,2] >= H1X[k,2] && HIV[i,2] <= H1X[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the TAD
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, H1X[k,2])
        EndCoord <- c(EndCoord, H1X[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_H1X <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_H1X <- HIV_H1X[order(HIV_H1X$Chromosome,HIV_H1X$StartCoord),]
head(HIV_H1X)
save(HIV_H1X, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_H1X.RData")

#############################
###     HIV_LIKE_H1.2     ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(LIKE_H1.2))) {
    if (HIV[i,1] == LIKE_H1.2[k,1] && HIV[i,2] >= LIKE_H1.2[k,2] && HIV[i,2] <= LIKE_H1.2[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the TAD
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, LIKE_H1.2[k,2])
        EndCoord <- c(EndCoord, LIKE_H1.2[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_LIKE_H1.2 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_LIKE_H1.2 <- HIV_LIKE_H1.2[order(HIV_LIKE_H1.2$Chromosome,HIV_LIKE_H1.2$StartCoord),]
head(HIV_LIKE_H1.2)
save(HIV_LIKE_H1.2, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_LIKE_H1.2.RData")

#############################
###     HIV_LIKE_H1X      ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV))) {
  for (k in seq(nrow(LIKE_H1X))) {
    if (HIV[i,1] == LIKE_H1X[k,1] && HIV[i,2] >= LIKE_H1X[k,2] && HIV[i,2] <= LIKE_H1X[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the TAD
      if (i > 1 && HIV[i,1] == HIV[i-1,1] && HIV[i,2] == HIV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HIV[i,1])
        InsertPos <- c(InsertPos,HIV[i,2])
        StartCoord <- c(StartCoord, LIKE_H1X[k,2])
        EndCoord <- c(EndCoord, LIKE_H1X[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HIV_LIKE_H1X <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HIV_LIKE_H1X <- HIV_LIKE_H1X[order(HIV_LIKE_H1X$Chromosome,HIV_LIKE_H1X$StartCoord),]
head(HIV_LIKE_H1X)
save(HIV_LIKE_H1X, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV_LIKE_H1X.RData")

#############################
###       HTLV_Gneg1      ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(Gneg1))) {
    if (HTLV[i,1] == Gneg1[k,1] && HTLV[i,2] >= Gneg1[k,2] && HTLV[i,2] <= Gneg1[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, Gneg1[k,2])
        EndCoord <- c(EndCoord, Gneg1[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_Gneg1 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_Gneg1 <- HTLV_Gneg1[order(HTLV_Gneg1$Chromosome,HTLV_Gneg1$StartCoord),]
save(HTLV_Gneg1, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_Gneg1.RData")
head(HTLV_Gneg1)

#############################
###       HTLV_Gneg2      ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(Gneg2))) {
    if (HTLV[i,1] == Gneg2[k,1] && HTLV[i,2] >= Gneg2[k,2] && HTLV[i,2] <= Gneg2[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, Gneg2[k,2])
        EndCoord <- c(EndCoord, Gneg2[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_Gneg2 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_Gneg2 <- HTLV_Gneg2[order(HTLV_Gneg2$Chromosome,HTLV_Gneg2$StartCoord),]
head(HTLV_Gneg2)
save(HTLV_Gneg2, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_Gneg2.RData")



#############################
###       HTLV_Gneg3     ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(Gneg3))) {
    if (HTLV[i,1] == Gneg3[k,1] && HTLV[i,2] >= Gneg3[k,2] && HTLV[i,2] <= Gneg3[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, Gneg3[k,2])
        EndCoord <- c(EndCoord, Gneg3[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_Gneg3 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_Gneg3 <- HTLV_Gneg3[order(HTLV_Gneg3$Chromosome,HTLV_Gneg3$StartCoord),]
head(HTLV_Gneg3)
save(HTLV_Gneg3, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_Gneg3.RData")


#############################
###       HTLV_Gneg4      ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(Gneg4))) {
    if (HTLV[i,1] == Gneg4[k,1] && HTLV[i,2] >= Gneg4[k,2] && HTLV[i,2] <= Gneg4[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, Gneg4[k,2])
        EndCoord <- c(EndCoord, Gneg4[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_Gneg4 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_Gneg4 <- HTLV_Gneg4[order(HTLV_Gneg4$Chromosome,HTLV_Gneg4$StartCoord),]
head(HTLV_Gneg4)
save(HTLV_Gneg4, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_Gneg4.RData")

#############################
###      HTLV_Gpos25      ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(Gpos25))) {
    if (HTLV[i,1] == Gpos25[k,1] && HTLV[i,2] >= Gpos25[k,2] && HTLV[i,2] <= Gpos25[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, Gpos25[k,2])
        EndCoord <- c(EndCoord, Gpos25[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_Gpos25 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_Gpos25 <- HTLV_Gpos25[order(HTLV_Gpos25$Chromosome,HTLV_Gpos25$StartCoord),]
head(HTLV_Gpos25)
save(HTLV_Gpos25, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_Gpos25.RData")


#############################
###      HTLV_Gpos50      ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(Gpos50))) {
    if (HTLV[i,1] == Gpos50[k,1] && HTLV[i,2] >= Gpos50[k,2] && HTLV[i,2] <= Gpos50[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, Gpos50[k,2])
        EndCoord <- c(EndCoord, Gpos50[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_Gpos50 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_Gpos50 <- HTLV_Gpos50[order(HTLV_Gpos50$Chromosome,HTLV_Gpos50$StartCoord),]
head(HTLV_Gpos50)
save(HTLV_Gpos50, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_Gpos50.RData")

#############################
###      HTLV_Gpos75      ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(Gpos75))) {
    if (HTLV[i,1] == Gpos75[k,1] && HTLV[i,2] >= Gpos75[k,2] && HTLV[i,2] <= Gpos75[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, Gpos75[k,2])
        EndCoord <- c(EndCoord, Gpos75[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_Gpos75 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_Gpos75 <- HTLV_Gpos75[order(HTLV_Gpos75$Chromosome,HTLV_Gpos75$StartCoord),]
head(HTLV_Gpos75)
save(HTLV_Gpos75, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_Gpos75.RData")


#############################
###      HTLV_Gpos100     ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(Gpos100))) {
    if (HTLV[i,1] == Gpos100[k,1] && HTLV[i,2] >= Gpos100[k,2] && HTLV[i,2] <= Gpos100[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the G Band
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, Gpos100[k,2])
        EndCoord <- c(EndCoord, Gpos100[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_Gpos100 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_Gpos100 <- HTLV_Gpos100[order(HTLV_Gpos100$Chromosome,HTLV_Gpos100$StartCoord),]
head(HTLV_Gpos100)
save(HTLV_Gpos100, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_Gpos100.RData")


#############################
###         HTLV_A        ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(A))) {
    if (HTLV[i,1] == A[k,1] && HTLV[i,2] >= A[k,2] && HTLV[i,2] <= A[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the A compartment
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, A[k,2])
        EndCoord <- c(EndCoord, A[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_A <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_A <- HTLV_A[order(HTLV_A$Chromosome,HTLV_A$StartCoord),]
head(HTLV_A)
save(HTLV_A, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_A.RData")

#############################
###         HTLV_B        ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(B))) {
    if (HTLV[i,1] == B[k,1] && HTLV[i,2] >= B[k,2] && HTLV[i,2] <= B[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the B compartment
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, B[k,2])
        EndCoord <- c(EndCoord, B[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_B <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_B <- HTLV_B[order(HTLV_B$Chromosome,HTLV_B$StartCoord),]
head(HTLV_B)
save(HTLV_B, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_B.RData")


#############################
###       HTLV_H1.2       ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(H1.2))) {
    if (HTLV[i,1] == H1.2[k,1] && HTLV[i,2] >= H1.2[k,2] && HTLV[i,2] <= H1.2[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the TAD
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, H1.2[k,2])
        EndCoord <- c(EndCoord, H1.2[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_H1.2 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_H1.2 <- HTLV_H1.2[order(HTLV_H1.2$Chromosome,HTLV_H1.2$StartCoord),]
head(HTLV_H1.2)
save(HTLV_H1.2, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_H1.2.RData")

#############################
###       HTLV_H1X        ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(H1X))) {
    if (HTLV[i,1] == H1X[k,1] && HTLV[i,2] >= H1X[k,2] && HTLV[i,2] <= H1X[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the TAD
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, H1X[k,2])
        EndCoord <- c(EndCoord, H1X[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_H1X <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_H1X <- HTLV_H1X[order(HTLV_H1X$Chromosome,HTLV_H1X$StartCoord),]
head(HTLV_H1X)
save(HTLV_H1X, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_H1X.RData")

#############################
###     HTLV_LIKE_H1.2    ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(LIKE_H1.2))) {
    if (HTLV[i,1] == LIKE_H1.2[k,1] && HTLV[i,2] >= LIKE_H1.2[k,2] && HTLV[i,2] <= LIKE_H1.2[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the TAD
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, LIKE_H1.2[k,2])
        EndCoord <- c(EndCoord, LIKE_H1.2[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_LIKE_H1.2 <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_LIKE_H1.2 <- HTLV_LIKE_H1.2[order(HTLV_LIKE_H1.2$Chromosome,HTLV_LIKE_H1.2$StartCoord),]
head(HTLV_LIKE_H1.2)
save(HTLV_LIKE_H1.2, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_LIKE_H1.2.RData")

#############################
###     HTLV_LIKE_H1X     ###
#############################

Chromosome <- as.numeric(as.character(c()))
InsertPos <- as.numeric(as.character(c()))
StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HTLV))) {
  for (k in seq(nrow(LIKE_H1X))) {
    if (HTLV[i,1] == LIKE_H1X[k,1] && HTLV[i,2] >= LIKE_H1X[k,2] && HTLV[i,2] <= LIKE_H1X[k,3]){ # if the chromosome is the same and the insert position falls inside the start and end coordinate of the TAD
      if (i > 1 && HTLV[i,1] == HTLV[i-1,1] && HTLV[i,2] == HTLV[i-1,2]){
        last <- tail(Count,n=1)
        Count <- head(Count,-1)
        counter <- last + 1
        Count <- c(Count,counter)
      }
      else {
        Chromosome <- c(Chromosome,HTLV[i,1])
        InsertPos <- c(InsertPos,HTLV[i,2])
        StartCoord <- c(StartCoord, LIKE_H1X[k,2])
        EndCoord <- c(EndCoord, LIKE_H1X[k,3])
        counter <- counter + 1
        Count <- c(Count,counter)
      }
    }
    else {
      counter <- 0
    }
  }
}

HTLV_LIKE_H1X <- data.frame(cbind(Chromosome,InsertPos,StartCoord,EndCoord,Count))
HTLV_LIKE_H1X <- HTLV_LIKE_H1X[order(HTLV_LIKE_H1X$Chromosome,HTLV_LIKE_H1X$StartCoord),]
head(HTLV_LIKE_H1X)
save(HTLV_LIKE_H1X, file="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HTLV_LIKE_H1X.RData")
