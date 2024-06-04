#"############ Starting bed creator script at $(date +'%H:%M:%S')... ##############"
##Muestra humanas normales no A431
#cd /home/vant/ATAC_seq/out
#export WD=$(pwd)
#for sample in $(ls | egrep _ )
#do
#    mkdir -p /home/vant/ATAC_seq/peak_calling/$sample
#    cd /home/vant/ATAC_seq/out/$sample/peak_calling
#    cut -f1-4 ATAC_peaks.broadPeak > $sample"_peaks.bed"
#    mv $sample*.bed /home/vant/ATAC_seq/peak_calling/$sample/
#done 
#
#cd /home/vant/ATAC_seq/peak_calling
#for sample in $(ls | egrep _O)
#do 
#    cat $sample/*.bed >> merge.bed
#done

#sort -k1,1 -k2,2n merge.bed > merge_sorted.bed

#bedtools merge -i merge_sorted.bed > merge_final.bed   

#"############ Finish bed files creation $(date +'%H:%M:%S')... ##############"

##Muestras humanas A431
#cd /home/vant/ATAC_seq/A431
#for sample in $(ls | egrep _ )
#do
#    mkdir -p /home/vant/ATAC_seq/peak_calling/A431/$sample
#    cd /home/vant/ATAC_seq/A431/$sample/peak_calling
#    cut -f1-4 ATAC_peaks.broadPeak > $sample"_peaks.bed"
#    mv $sample*.bed /home/vant/ATAC_seq/peak_calling/A431/$sample/
#done 
#
#cd /home/vant/ATAC_seq/peak_calling/A431
#for sample in $(ls | egrep _)
#do 
#    cat $sample/*.bed >> merge_A431.bed
#done

#sort -k1,1 -k2,2n merge_A431.bed > merge_sorted_A431.bed

#bedtools merge -i merge_sorted_A431.bed > merge_final_A431.bed   

##Muestras de ratón
cd /home/vant/ATAC_seq/mouse_ATAC
for sample in $(ls | egrep _ )
do
    mkdir -p /home/vant/ATAC_seq/peak_calling/mouse_data/$sample
    cd /home/vant/ATAC_seq/mouse_ATAC/$sample/peak_calling
    cut -f1-4 ATAC_peaks.broadPeak > $sample"_peaks.bed"
    mv $sample*.bed /home/vant/ATAC_seq/peak_calling/mouse_data/$sample/
done 

cd /home/vant/ATAC_seq/peak_calling/mouse_data
for sample in $(ls | egrep s1)
do 
    cat $sample/*.bed >> merge_mouse.bed
done

sort -k1,1 -k2,2n merge_mouse.bed > merge_sorted_mouse.bed

bedtools merge -i merge_sorted_mouse.bed > merge_final_mouse.bed   