#Script to create exonic, intronic, and intergenic beds.
#Adapted from http://crazyhottommy.blogspot.com/2013/05/find-exons-introns-and-intergenic.html

GTF_FILE=$1
EXON_BED=$2
INTRON_BED=$3

#Create Exon Bed
cat "${GTF_FILE}" | awk 'BEGIN{FS="\t";OFS="\t";} $3=="exon" {print $1,$4,$5,$9}' > "${EXON_BED}"

sortBed -i "${EXON_BED}" > tmp.bed

mv -f tmp.bed "${EXON_BED}"

#Merge Exons Together

mergeBed -i "${EXON_BED}" -c 1,4 -o count,collapse > tmp_merged_exon.bed
mv -f tmp_merged_exon.bed "${EXON_BED}"

#Create Intron Bed
cat "${GTF_FILE}" | awk 'BEGIN{FS="\t";OFS="\t";} $3=="gene" {print $1,$4,$5,$9}' > tmp.bed

sortBed -i tmp.bed > tmp_gene.bed

subtractBed -a tmp_gene.bed -b "${EXON_BED}" > "${INTRON_BED}"
rm tmp.bed
rm tmp_gene.bed
