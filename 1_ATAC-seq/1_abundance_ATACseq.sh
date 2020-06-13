# Mapping the abundance of ATAC-seq in all genomic structures

##################
#    GBands      #
##################

bedtools map -a Gneg1.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_Gneg1.bed
bedtools map -a Gneg2.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_Gneg2.bed
bedtools map -a Gneg3.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_Gneg3.bed
bedtools map -a Gneg4.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_Gneg4.bed
bedtools map -a Gpos25.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_Gpos25.bed
bedtools map -a Gpos50.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_Gpos50.bed
bedtools map -a Gpos75.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_Gpos75.bed
bedtools map -a Gpos100.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_Gpos100.bed


#################################
#       A/B compartments        #
#################################

bedtools sort -i Acompartment.bed > sorted_A.bed
bedtools map -a sorted_A.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_A.bed

bedtools sort -i Bcompartment.bed > sorted_B.bed
bedtools map -a sorted_B.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_B.bed


##################
#      TADs      #
##################

bedtools map -a H1.2_enriched.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_H1.2.bed
bedtools map -a H1X_enriched.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_H1X.bed
bedtools map -a LIKE_H1.2_enriched.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_LIKE_H1.2.bed
bedtools map -a LIKE_H1X_enriched.bed -b ATACSeq_WT_100kb.bed -c 4 -o mean > ATACSeq_LIKE_H1X.bed



