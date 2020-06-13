# Mapping the abundance of histones H1.2 and H1X in all genomic structures

# Gneg1

bedtools map -a Gneg1.bed -b H1.2endo.bed -c 5 -o mean > H1.2_Gneg1.bed
bedtools map -a Gneg1.bed -b H1Xendo.bed -c 5 -o mean > H1X_Gneg1.bed


# Gneg2

bedtools map -a Gneg2.bed -b H1.2endo.bed -c 5 -o mean > H1.2_Gneg2.bed
bedtools map -a Gneg2.bed -b H1Xendo.bed -c 5 -o mean > H1X_Gneg2.bed


# Gneg3

bedtools map -a Gneg3.bed -b H1.2endo.bed -c 5 -o mean > H1.2_Gneg3.bed
bedtools map -a Gneg3.bed -b H1Xendo.bed -c 5 -o mean > H1X_Gneg3.bed


# Gneg4

bedtools map -a Gneg4.bed -b H1.2endo.bed -c 5 -o mean > H1.2_Gneg4.bed
bedtools map -a Gneg4.bed -b H1Xendo.bed -c 5 -o mean > H1X_Gneg4.bed


# Gpos25

bedtools map -a Gpos25.bed -b H1.2endo.bed -c 5 -o mean > H1.2_Gpos25.bed
bedtools map -a Gpos25.bed -b H1Xendo.bed -c 5 -o mean > H1X_Gpos25.bed


# Gpos50

bedtools map -a Gpos50.bed -b H1.2endo.bed -c 5 -o mean > H1.2_Gpos50.bed
bedtools map -a Gpos50.bed -b H1Xendo.bed -c 5 -o mean > H1X_Gpos50.bed


# Gpos75

bedtools map -a Gpos75.bed -b H1.2endo.bed -c 5 -o mean > H1.2_Gpos75.bed
bedtools map -a Gpos75.bed -b H1Xendo.bed -c 5 -o mean > H1X_Gpos75.bed


# Gpos100

bedtools map -a Gpos100.bed -b H1.2endo.bed -c 5 -o mean > H1.2_Gpos100.bed
bedtools map -a Gpos100.bed -b H1Xendo.bed -c 5 -o mean > H1X_Gpos100.bed


# A compartment

bedtools sort -i Acompartment.bed > sorted_A.bed
bedtools map -a sorted_A.bed -b H1.2endo.bed -c 5 -o mean > H1.2_A.bed
bedtools map -a sorted_A.bed -b H1Xendo.bed -c 5 -o mean > H1X_A.bed


# B compartment

bedtools sort -i Bcompartment.bed > sorted_B.bed
bedtools map -a sorted_B.bed -b H1.2endo.bed -c 5 -o mean > H1.2_B.bed
bedtools map -a sorted_B.bed -b H1Xendo.bed -c 5 -o mean > H1X_B.bed


# H1.2 enriched

bedtools map -a H1.2_enriched.bed -b H1.2endo.bed -c 5 -o mean > H1.2_H1.2.bed
bedtools map -a H1.2_enriched.bed -b H1Xendo.bed -c 5 -o mean > H1X_H1.2.bed


# H1X enriched

bedtools map -a H1X_enriched.bed -b H1.2endo.bed -c 5 -o mean > H1.2_H1X.bed
bedtools map -a H1X_enriched.bed -b H1Xendo.bed -c 5 -o mean > H1X_H1X.bed


# LIKE_H1.2 enriched

bedtools map -a LIKE_H1.2_enriched.bed -b H1.2endo.bed -c 5 -o mean > H1.2_LIKE_H1.2.bed
bedtools map -a LIKE_H1.2_enriched.bed -b H1Xendo.bed -c 5 -o mean > H1X_LIKE_H1.2.bed


# LIKE_H1X enriched

bedtools map -a LIKE_H1X_enriched.bed -b H1.2endo.bed -c 5 -o mean > H1.2_LIKE_H1X.bed
bedtools map -a LIKE_H1X_enriched.bed -b H1Xendo.bed -c 5 -o mean > H1X_LIKE_H1X.bed
