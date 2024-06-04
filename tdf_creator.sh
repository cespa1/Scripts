cd /home/vant/ATAC_seq/mouse_ATAC
for muestra in $(ls | grep s1)
do
    igvtools count -z 5 -w 25 -e 250 $muestra/$muestra*.bam $muestra/$muestra".tdf" /home/vant/Descargas/mm10.fa
done
