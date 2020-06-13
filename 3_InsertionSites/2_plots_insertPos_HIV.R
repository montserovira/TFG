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


######################
###   HIV GBands   ###
######################

load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_Gneg1.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_Gneg2.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_Gneg3.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_Gneg4.RData")

load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_Gpos25.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_Gpos50.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_Gpos75.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_Gpos100.RData")


#############################
###       HIV_Gneg1       ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_Gneg1))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_Gneg1[i,3])
    EndCoord <- c(EndCoord, HIV_Gneg1[i,4])
    counter <- counter + HIV_Gneg1[i,5]
  }
  else if (i > 1 && HIV_Gneg1[i,3] == HIV_Gneg1[i-1,3]){
    counter <- counter + HIV_Gneg1[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_Gneg1[i,3])
    EndCoord <- c(EndCoord, HIV_Gneg1[i,4])
    counter <- counter + HIV_Gneg1[i,5]
  }
  if (i == nrow(HIV_Gneg1)){
    Count <- c(Count,counter)
  }
}

new_HIV_Gneg1 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_Gneg1$Ratio <- ((new_HIV_Gneg1$Count)/(new_HIV_Gneg1$EndCoord - new_HIV_Gneg1$StartCoord))
head(new_HIV_Gneg1)


#############################
###       HIV_Gneg2       ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_Gneg2))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_Gneg2[i,3])
    EndCoord <- c(EndCoord, HIV_Gneg2[i,4])
    counter <- counter + HIV_Gneg2[i,5]
  }
  else if (i > 1 && HIV_Gneg2[i,3] == HIV_Gneg2[i-1,3]){
    counter <- counter + HIV_Gneg2[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_Gneg2[i,3])
    EndCoord <- c(EndCoord, HIV_Gneg2[i,4])
    counter <- counter + HIV_Gneg2[i,5]
  }
  if (i == nrow(HIV_Gneg2)){
    Count <- c(Count,counter)
  }
}

new_HIV_Gneg2 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_Gneg2$Ratio <- ((new_HIV_Gneg2$Count)/(new_HIV_Gneg2$EndCoord - new_HIV_Gneg2$StartCoord))
head(new_HIV_Gneg2)


#############################
###       HIV_Gneg3       ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_Gneg3))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_Gneg3[i,3])
    EndCoord <- c(EndCoord, HIV_Gneg3[i,4])
    counter <- counter + HIV_Gneg3[i,5]
  }
  else if (i > 1 && HIV_Gneg3[i,3] == HIV_Gneg3[i-1,3]){
    counter <- counter + HIV_Gneg3[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_Gneg3[i,3])
    EndCoord <- c(EndCoord, HIV_Gneg3[i,4])
    counter <- counter + HIV_Gneg3[i,5]
  }
  if (i == nrow(HIV_Gneg3)){
    Count <- c(Count,counter)
  }
}

new_HIV_Gneg3 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_Gneg3$Ratio <- ((new_HIV_Gneg3$Count)/(new_HIV_Gneg3$EndCoord - new_HIV_Gneg3$StartCoord))
head(new_HIV_Gneg3)


#############################
###       HIV_Gneg4       ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_Gneg4))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_Gneg4[i,3])
    EndCoord <- c(EndCoord, HIV_Gneg4[i,4])
    counter <- counter + HIV_Gneg4[i,5]
  }
  else if (i > 1 && HIV_Gneg4[i,3] == HIV_Gneg4[i-1,3]){
    counter <- counter + HIV_Gneg4[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_Gneg4[i,3])
    EndCoord <- c(EndCoord, HIV_Gneg4[i,4])
    counter <- counter + HIV_Gneg4[i,5]
  }
  if (i == nrow(HIV_Gneg4)){
    Count <- c(Count,counter)
  }
}

new_HIV_Gneg4 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_Gneg4$Ratio <- ((new_HIV_Gneg4$Count)/(new_HIV_Gneg4$EndCoord - new_HIV_Gneg4$StartCoord))
head(new_HIV_Gneg4)


#############################
###       HIV_Gpos25      ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_Gpos25))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_Gpos25[i,3])
    EndCoord <- c(EndCoord, HIV_Gpos25[i,4])
    counter <- counter + HIV_Gpos25[i,5]
  }
  else if (i > 1 && HIV_Gpos25[i,3] == HIV_Gpos25[i-1,3]){
    counter <- counter + HIV_Gpos25[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_Gpos25[i,3])
    EndCoord <- c(EndCoord, HIV_Gpos25[i,4])
    counter <- counter + HIV_Gpos25[i,5]
  }
  if (i == nrow(HIV_Gpos25)){
    Count <- c(Count,counter)
  }
}

new_HIV_Gpos25 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_Gpos25$Ratio <- ((new_HIV_Gpos25$Count)/(new_HIV_Gpos25$EndCoord - new_HIV_Gpos25$StartCoord))
head(new_HIV_Gpos25)


#############################
###       HIV_Gpos50      ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_Gpos50))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_Gpos50[i,3])
    EndCoord <- c(EndCoord, HIV_Gpos50[i,4])
    counter <- counter + HIV_Gpos50[i,5]
  }
  else if (i > 1 && HIV_Gpos50[i,3] == HIV_Gpos50[i-1,3]){
    counter <- counter + HIV_Gpos50[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_Gpos50[i,3])
    EndCoord <- c(EndCoord, HIV_Gpos50[i,4])
    counter <- counter + HIV_Gpos50[i,5]
  }
  if (i == nrow(HIV_Gpos50)){
    Count <- c(Count,counter)
  }
}

new_HIV_Gpos50 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_Gpos50$Ratio <- ((new_HIV_Gpos50$Count)/(new_HIV_Gpos50$EndCoord - new_HIV_Gpos50$StartCoord))
head(new_HIV_Gpos50)


#############################
###       HIV_Gpos75      ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_Gpos75))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_Gpos75[i,3])
    EndCoord <- c(EndCoord, HIV_Gpos75[i,4])
    counter <- counter + HIV_Gpos75[i,5]
  }
  else if (i > 1 && HIV_Gpos75[i,3] == HIV_Gpos75[i-1,3]){
    counter <- counter + HIV_Gpos75[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_Gpos75[i,3])
    EndCoord <- c(EndCoord, HIV_Gpos75[i,4])
    counter <- counter + HIV_Gpos75[i,5]
  }
  if (i == nrow(HIV_Gpos75)){
    Count <- c(Count,counter)
  }
}

new_HIV_Gpos75 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_Gpos75$Ratio <- ((new_HIV_Gpos75$Count)/(new_HIV_Gpos75$EndCoord - new_HIV_Gpos75$StartCoord))
head(new_HIV_Gpos75)


#############################
###      HIV_Gpos100      ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_Gpos100))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_Gpos100[i,3])
    EndCoord <- c(EndCoord, HIV_Gpos100[i,4])
    counter <- counter + HIV_Gpos100[i,5]
  }
  else if (i > 1 && HIV_Gpos100[i,3] == HIV_Gpos100[i-1,3]){
    counter <- counter + HIV_Gpos100[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_Gpos100[i,3])
    EndCoord <- c(EndCoord, HIV_Gpos100[i,4])
    counter <- counter + HIV_Gpos100[i,5]
  }
  if (i == nrow(HIV_Gpos100)){
    Count <- c(Count,counter)
  }
}

new_HIV_Gpos100 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_Gpos100$Ratio <- ((new_HIV_Gpos100$Count)/(new_HIV_Gpos100$EndCoord - new_HIV_Gpos100$StartCoord))
head(new_HIV_Gpos100)


GBands_Colors <- c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc","Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3")

boxplot(new_HIV_Gneg1$Ratio,new_HIV_Gneg2$Ratio,new_HIV_Gneg3$Ratio,new_HIV_Gneg4$Ratio,
        new_HIV_Gpos25$Ratio,new_HIV_Gpos50$Ratio,new_HIV_Gpos75$Ratio,new_HIV_Gpos100$Ratio,
        col = GBands_Colors, outline = F, main = "HIV GBands Insert Positions (normalyzed by length)", 
        ylab="Ratio",  at=c(1,3,5,7,9,11,13,15), xaxt="n")

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


######################
###     HIV AB     ###
######################

load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_A.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_B.RData")


#############################
###         HIV_A         ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_A))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_A[i,3])
    EndCoord <- c(EndCoord, HIV_A[i,4])
    counter <- counter + HIV_A[i,5]
  }
  else if (i > 1 && HIV_A[i,3] == HIV_A[i-1,3]){
    counter <- counter + HIV_A[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_A[i,3])
    EndCoord <- c(EndCoord, HIV_A[i,4])
    counter <- counter + HIV_A[i,5]
  }
  if (i == nrow(HIV_A)){
    Count <- c(Count,counter)
  }
}

new_HIV_A <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_A$Ratio <- ((new_HIV_A$Count)/(new_HIV_A$EndCoord - new_HIV_A$StartCoord))
head(new_HIV_A)


#############################
###         HIV_B         ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_B))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_B[i,3])
    EndCoord <- c(EndCoord, HIV_B[i,4])
    counter <- counter + HIV_B[i,5]
  }
  else if (i > 1 && HIV_B[i,3] == HIV_B[i-1,3]){
    counter <- counter + HIV_B[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_B[i,3])
    EndCoord <- c(EndCoord, HIV_B[i,4])
    counter <- counter + HIV_B[i,5]
  }
  if (i == nrow(HIV_B)){
    Count <- c(Count,counter)
  }
}

new_HIV_B <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_B$Ratio <- ((new_HIV_B$Count)/(new_HIV_B$EndCoord - new_HIV_B$StartCoord))
head(new_HIV_B)


AB_Colors <- c("A_compartment"="sandybrown","B_compartment"="lightblue2")

boxplot(new_HIV_A$Ratio,new_HIV_B$Ratio,
        col = AB_Colors, outline = F, main = "HIV A/B compartments Insert Positions (normalyzed by length)", 
        ylab="Ratio",  at=c(1,3), xaxt="n")

abline(v=2, col = "grey", lwd = 1)

title(sub="A compartment", adj=0.2, line=0.8)
title(sub="B compartment", adj=0.85, line=0.8)


######################
###    HIV TADs    ###
######################

load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_H1.2.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_H1X.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_LIKE_H1.2.RData")
load("/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/HIV/HIV_LIKE_H1X.RData")


#############################
###       HIV_H1.2        ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_H1.2))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_H1.2[i,3])
    EndCoord <- c(EndCoord, HIV_H1.2[i,4])
    counter <- counter + HIV_H1.2[i,5]
  }
  else if (i > 1 && HIV_H1.2[i,3] == HIV_H1.2[i-1,3]){
    counter <- counter + HIV_H1.2[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_H1.2[i,3])
    EndCoord <- c(EndCoord, HIV_H1.2[i,4])
    counter <- counter + HIV_H1.2[i,5]
  }
  if (i == nrow(HIV_H1.2)){
    Count <- c(Count,counter)
  }
}

new_HIV_H1.2 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_H1.2$Ratio <- ((new_HIV_H1.2$Count)/(new_HIV_H1.2$EndCoord - new_HIV_H1.2$StartCoord))
head(new_HIV_H1.2)


#############################
###        HIV_H1X        ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_H1X))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_H1X[i,3])
    EndCoord <- c(EndCoord, HIV_H1X[i,4])
    counter <- counter + HIV_H1X[i,5]
  }
  else if (i > 1 && HIV_H1X[i,3] == HIV_H1X[i-1,3]){
    counter <- counter + HIV_H1X[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_H1X[i,3])
    EndCoord <- c(EndCoord, HIV_H1X[i,4])
    counter <- counter + HIV_H1X[i,5]
  }
  if (i == nrow(HIV_H1X)){
    Count <- c(Count,counter)
  }
}

new_HIV_H1X <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_H1X$Ratio <- ((new_HIV_H1X$Count)/(new_HIV_H1X$EndCoord - new_HIV_H1X$StartCoord))
head(new_HIV_H1X)


#############################
###   HIV_LIKE_H1.2      ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_LIKE_H1.2))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_LIKE_H1.2[i,3])
    EndCoord <- c(EndCoord, HIV_LIKE_H1.2[i,4])
    counter <- counter + HIV_LIKE_H1.2[i,5]
  }
  else if (i > 1 && HIV_LIKE_H1.2[i,3] == HIV_LIKE_H1.2[i-1,3]){
    counter <- counter + HIV_LIKE_H1.2[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_LIKE_H1.2[i,3])
    EndCoord <- c(EndCoord, HIV_LIKE_H1.2[i,4])
    counter <- counter + HIV_LIKE_H1.2[i,5]
  }
  if (i == nrow(HIV_LIKE_H1.2)){
    Count <- c(Count,counter)
  }
}

new_HIV_LIKE_H1.2 <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_LIKE_H1.2$Ratio <- ((new_HIV_LIKE_H1.2$Count)/(new_HIV_LIKE_H1.2$EndCoord - new_HIV_LIKE_H1.2$StartCoord))
head(new_HIV_LIKE_H1.2)


#############################
###     HIV_LIKE_H1X      ###
#############################

StartCoord <- as.numeric(as.character(c()))
EndCoord <- as.numeric(as.character(c()))
counter <- 0
Count <- as.numeric(as.character(c()))

for (i in seq(nrow(HIV_LIKE_H1X))){
  if (i == 1){
    StartCoord <- c(StartCoord, HIV_LIKE_H1X[i,3])
    EndCoord <- c(EndCoord, HIV_LIKE_H1X[i,4])
    counter <- counter + HIV_LIKE_H1X[i,5]
  }
  else if (i > 1 && HIV_LIKE_H1X[i,3] == HIV_LIKE_H1X[i-1,3]){
    counter <- counter + HIV_LIKE_H1X[i,5]
  }
  else {
    Count <- c(Count,counter)
    counter <- 0
    StartCoord <- c(StartCoord, HIV_LIKE_H1X[i,3])
    EndCoord <- c(EndCoord, HIV_LIKE_H1X[i,4])
    counter <- counter + HIV_LIKE_H1X[i,5]
  }
  if (i == nrow(HIV_LIKE_H1X)){
    Count <- c(Count,counter)
  }
}

new_HIV_LIKE_H1X <- data.frame(cbind(StartCoord,EndCoord,Count))
new_HIV_LIKE_H1X$Ratio <- ((new_HIV_LIKE_H1X$Count)/(new_HIV_LIKE_H1X$EndCoord - new_HIV_LIKE_H1X$StartCoord))
head(new_HIV_LIKE_H1X)


TADs_Colors <- c("H1X_enriched"="sandybrown","LIKE_H1X_enriched"="coral2","H1.2_enriched"="lightblue2","LIKE_H1.2_enriched"="cyan3")

boxplot(new_HIV_H1X$Ratio,new_HIV_LIKE_H1X$Ratio,new_HIV_LIKE_H1.2$Ratio,new_HIV_H1.2$Ratio,
        col = TADs_Colors, outline = F, main = "HIV TADs Insert Positions (normalyzed by length)", 
        ylab="Ratio",  at=c(1,3,5,7), xaxt="n")

abline(v=2, col = "grey", lwd = 1)
abline(v=4, col = "grey", lwd = 1)
abline(v=6, col = "grey", lwd = 1)

title(sub="H1X", adj=0.1, line=0.8)
title(sub="LIKE_H1X", adj=0.36, line=0.8)
title(sub="LIKE_H1.2", adj=0.65, line=0.8)
title(sub="H1.2", adj=0.91, line=0.8)
