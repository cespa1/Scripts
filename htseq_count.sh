#cd /home/vant/ATAC_seq/peak_calling
#awk '{print $1"\tunknown\texon\t"$2"\t"$3"\t.\t+\t0\tgene_id \42"$1":"$2"-"$3"\42;\ttranscript_id \42"$1":"$2"-"$3"\42;\tgene_name \42"$1":"$2"-"$3"\42;"}' merge_final.bed  > all_peaks_merged.gtf
#awk '{print $1"\tunknown\texon\t"$2"\t"$3"\t.\t+\t0\tgene_id \42"$1":"$2"-"$3"\42;\ttranscript_id \42"$1":"$2"-"$3"\42;\tgene_name \42"$1":"$2"-"$3"\42;"}' A431/merge_final_A431.bed  > A431/all_peaks_merged_A431.gtf
#cd /home/vant/ATAC_seq/peak_calling/mouse_data
#awk '{print $1"\tunknown\texon\t"$2"\t"$3"\t.\t+\t0\tgene_id \42"$1":"$2"-"$3"\42;\ttranscript_id \42"$1":"$2"-"$3"\42;\tgene_name \42"$1":"$2"-"$3"\42;"}' merge_final_mouse.bed > all_peaks_merged_mouse.gtf

#cd /home/vant/ATAC_seq
#for muestra in $(ls out/ | egrep _O)
#do
#    htseq-count --format=bam --type=exon --idattr=gene_name out/$muestra/*.bam peak_calling/all_peaks_merged.gtf  > peak_calling/$muestra/$muestra"_htseq_counts"
#done
#A431
#cd /home/vant/ATAC_seq
#for muestra in $(ls A431/ | egrep .bam | cut -d'_' -f1,2,3,4)
#do
#    htseq-count --format=bam --type=exon --idattr=gene_name A431/$muestra*.bam A431/peak_calling/all_peaks_merged_A431.gtf  > A431/peak_calling/htseq/$muestra"_htseq_counts"
#done
#Mouse
cd /home/vant/ATAC_seq
for muestra in $(ls mouse_ATAC/ | egrep .bam | cut -d "." -f1) 
do
    htseq-count --format=bam --type=exon --idattr=gene_name mouse_ATAC/$muestra.bam mouse_ATAC/peak_calling/all_peaks_merged_mouse.gtf  > /home/vant/ATAC_seq/mouse_ATAC/peak_calling/htseq/$muestra"_htseq_counts"
done




