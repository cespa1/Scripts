#cd Mapping_test
#cd fastq
mkdir -p $1
cd $1
for muestra in $(ls | grep .bam | cut -d'_' -f1,2,3,4)
do
mkdir -p peak_calling
cd peak_calling 
touch $muestra"_ATAC_peaks.broadPeak"/$muestra"_peaks.bed"
cut -f1-4 $muestra"_ATAC_peaks.broadPeak"/*.broadPeak > $muestra"_ATAC_peaks.broadPeak"/$muestra"_peaks.bed"
cat $muestra"_ATAC_peaks.broadPeak"/$muestra"_peaks.bed" >> merge_peaks_$1.bed
cd ..
done
cd peak_calling
sort -k1,1 -k2,2n merge_peaks_$1.bed > merge_sorted_$1.bed

bedtools merge -i merge_sorted_$1.bed > merge_final_$1.bed   

awk '{print $1"\tunknown\texon\t"$2"\t"$3"\t.\t+\t0\tgene_id \42"$1":"$2"-"$3"\42;\ttranscript_id \42"$1":"$2"-"$3"\42;\tgene_name \42"$1":"$2"-"$3"\42;"}' merge_final_$1.bed    > all_peaks_merged_$1.gtf

cd ..

for muestra in $(ls | grep .bam | cut -d'.' -f1)
do
mkdir -p htseq_counts
    htseq-count --format=bam --type=exon --idattr=gene_name $muestra.bam peak_calling/all_peaks_merged_$1.gtf  > htseq_counts/$muestra"_htseq_counts"
done
