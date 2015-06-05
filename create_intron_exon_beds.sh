#Script to create exonic, intronic, and intergenic beds.
#Adapted from http://crazyhottommy.blogspot.com/2013/05/find-exons-introns-and-intergenic.html

GTF_FILE=$1
EXON_BED=$2
INTRON_BED=$3

#Create Exon Bed
cat $GTF_FILE | awk 'BEGIN{FS="\t";OFS="\t";} $3=="exon" {print $1,$4,$5,$9}' > $EXON_BED

sortBed -i $EXON_BED > temp.bed

mv -f temp.bed $EXON_BED

#Merge Exons Together

mergeBed -i $EXON_BED -c 1,4 -o count,collapse > tmp_merged_exon.bed


#Create Intron Bed
cat $GTF_FILE | awk 'BEGIN{FS="\t";OFS="\t";} $3=="gene" {print $1,$4,$5,$9}' > tmp.bed

sortBed -i tmp.bed > tmp_gene.bed

subtractBed -a tmp_gene.bed -b tmp_merged_exon.bed > $INTRON_BED
rm tmp_merged_exon.bed
rm tmp.bed
rm tmp_gene.bed
