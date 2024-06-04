##Peak calling con MACS2 HUMANO
#cd /home/vant/ATAC_seq/out
#for muestra in $(ls | grep _)
#do
#if [ -a $muestra/peak_calling/ATAC_peaks.broadPeak ]
#then
#  echo  "$muestra peak calling is alredy done going next sample"
#else
#   macs2 callpeak -t $muestra/*dedup.bam --outdir $muestra/peak_calling/ \
#        -f BAM -n ATAC -g hs -q 0.0001 --nomodel --shift -37 --extsize 73 --broad
#fi
#done

##PEAKCALLING RATON
cd /home/vant/ATAC_seq/mouse_ATAC
for muestra in $(ls | grep _)
do
if [ -a $muestra/peak_calling/ATAC_peaks.broadPeak ]
then
  echo  "$muestra peak calling is alredy done going next sample"
else
    macs2 callpeak -t $muestra/*dedup.bam --outdir $muestra/peak_calling/ \
        -f BAM -n ATAC -g mm -q 0.0001 --nomodel --shift -37 --extsize 73 --broad
fi
done

##PEAKCALLING A431 (Humano)  
#cd /home/vant/ATAC_seq/A431
#for muestra in $(ls | grep _)
#do
#if [ -a $muestra/peak_calling/ATAC_peaks.broadPeak ]
#then
#  echo  "$muestra peak calling is alredy done going next sample"
#else
#    macs2 callpeak -t $muestra/*dedup.bam --outdir $muestra/peak_calling/ \
##        -f BAM -n ATAC -g hs -q 0.0001 --nomodel --shift -37 --extsize 73 --broad
#fi
#done
