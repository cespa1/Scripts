
##Peak calling
cd Bam/
for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
        macs2 callpeak -t $muestra*dedup.bam --outdir ../peak_calling/$muestra"_ATAC_peaks.broadPeak" \
    -f BAM -n ATAC -g hs -q 0.0001 --nomodel --bw 1000 --shift -37 --extsize 73 --broad

done
cd ..

##Peakome

touch merge.bed
for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
cut -f1-4 $muestra"_ATAC_peaks.broadPeak"/ATAC_peaks.broadPeak > $muestra"_ATAC_peaks.broadPeak"/$muestra"_peaks.bed"
cat $muestra"_ATAC_peaks.broadPeak"/$muestra"_peaks.bed" >> merge.bed
done

sort -k1,1 -k2,2n merge.bed > sorted_merge.bed

bedtools merge -i sorted_merge.bed > merge_final.bed
rm merge.bed