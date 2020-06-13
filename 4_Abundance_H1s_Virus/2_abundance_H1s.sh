WorkingDir="/Users/montseroviraripoll/Documents/Bioinformatica/4🍟/TFG/Retrovirus/Abundance"
cd $WorkingDir

# convert WIG format to BED format

convert2bed -i wig < H1.2endo_r2.fastq.gz.bdg_InputSubtracted.wig > H1.2endo.bed
convert2bed -i wig < H1Xendo_r3.fastq.gz.bdg_InputSubtracted.wig > H1Xendo.bed

awk '{print "chr"$0}' HIV_bed.bed > HIV_bed.bed
awk '{print "chr"$0}' HTLV_bed.bed > HTLV_bed.bed

bedtools sort -i HIV_bed.bed > sorted_HIV.bed
bedtools sort -i HTLV_bed.bed > sorted_HTLV.bed

# Mapping the abundance of histones H1.2 and H1X in both retroviruses HIV and HTLV

bedtools map -a sorted_HIV.bed -b H1.2endo.bed -c 5 -o mean > H1.2_HIV.bed
bedtools map -a sorted_HTLV.bed -b H1.2endo.bed -c 5 -o mean > H1.2_HTLV.bed

bedtools map -a sorted_HIV.bed -b H1Xendo.bed -c 5 -o mean > H1X_HIV.bed
bedtools map -a sorted_HTLV.bed -b H1Xendo.bed -c 5 -o mean > H1X_HTLV.bed

