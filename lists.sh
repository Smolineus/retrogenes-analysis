### Creating 3 lists: intergenic, intronic and with 1 exon CDS ###

##INTERGENIC - not in the transcripts##

#Taking only 'transcript columns from adnotations for faster procesing
# źle zrobione bo MERGE: zcat hs1.ncbiRefSeq.bigZip.gtf.gz | awk '$3=="transcript" {print $1"\t"($4-1)"\t"$5}' | bedtools sort -i stdin | bedtools merge -i stdin > transcript_cordinates.bed

# awk 4-1 (bo pliki gtf są 1-based, a bed 0-based)
zcat hs1.ncbiRefSeq.bigZip.gtf.gz | awk '$3=="transcript" {print $1"\t"($4-1)"\t"$5"\t"$10}' > transcript_cordinates_all.bed
bedtools intersect -a human_retrocopies_merged.bed  -b transcript_cordinates_all.bed -v -wo > intergenic_retrogenes.bed

echo 'INTERGENIC DONE'

##INTRONIC - inside transcripts, but not in the exon##

zcat hs1.ncbiRefSeq.bigZip.gtf.gz | awk '$3 == "exon" {print $1"\t"$4-1"\t"$5"\t"$10}' > exons.bed

bedtools intersect -a human_retrocopies_merged.bed -b transcript_cordinates_all.bed -f 1.0 -wo > inside_transcripts.bed

bedtools intersect -a inside_transcripts.bed -b exons.bed -v -wo > intronic_retrogenes.bed

# filtrowanie - zostawiamy tylko interesujące nas linijki i zostawiamy gen 

cat intronic_retrogenes.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $(NF-1)}' intronic_retrogenes.bed > intronic_retrogenes_genes.bed

echo 'INTRONIC DONE'

## CDS = 1 exon ##

#single_cds.bed -> z mojego skryptu CDS_one_exon.py

bedtools intersect -a human_retrocopies_merged.bed -b single_cds.bed -wa -wb -f 1.0 > retrogenes_in_single_exon_cds.bed

echo 'CDS DONE'

###INTERSECT with encode long reads### - PÓŹNIEJ narazie same listy

#bedtools intersect -a encode4_long_liftover.bed -b intergenic_retrogenes.bed -wo -f 0.1 -F 0.1 > intergenic_encode.tsv

#bedtools intersect -a encode4_long_liftover.bed -b intronic_retrogenes.bed -wo -f 0.1 -F 0.1 > intronic_encode.tsv

#bedtools intersect -a encode4_long_liftover.bed -b retrogenes_in_single_exon_cds.bed -wo -f 0.1 -F 0.1 > single_exon_cds_encode.tsv

echo 'INTERSECTS DONE'
