##Bucle ATACorrect
echo "high or low?"
read high_low

until [ $high_low == "high" ] || [ $high_low ==  "low" ]
do 
    echo "high or low?"
    read nombre_muestra
done

if [ $high_low == "high" ]
then
for bam in $(ls | egrep bam | cut -d'_' -f1,2,3,4)
do

    TOBIAS ATACorrect --bam $bam*.bam --genome /home/vant/ATAC_seq/OSCC/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --peaks /home/vant/ATAC_seq/OSCC/prueba/beds/BATF3_high.bed --blacklist /home/vant/ATAC_seq/OSCC/footprinting/blacklist/hg38.blacklist.bed --outdir $bam/BATF3_high/ --cores 7

done

for bam in $(ls | egrep bam | cut -d'_' -f1,2,3,4)
do

    TOBIAS FootprintScores --signal $bam/BATF3_high/$bam*"_corrected.bw" --regions /home/vant/ATAC_seq/OSCC/prueba/beds/BATF3_high.bed \
    --output $bam/BATF3_high/$bam"_footprints_BATF3_high.bw" --cores 7

done

for bam in $(ls | egrep bam | cut -d'_' -f1,2,3,4)
do
    
    TOBIAS PlotAggregate --TFBS /home/vant/ATAC_seq/OSCC/prueba/beds/BATF3_high.bed --signals $bam/BATF3_high/$bam*"_corrected.bw" $bam/BATF3_high/$bam*"_footprints_BATF3_high.bw" --output $bam/BATF3_high/footprint.pdf --share_y both --plot_boundaries --signal-on-x --flank 250

done

else
for bam in $(ls | egrep bam | cut -d'_' -f1,2,3,4)
do

    TOBIAS ATACorrect --bam $bam*.bam --genome /home/vant/ATAC_seq/OSCC/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --peaks /home/vant/ATAC_seq/OSCC/prueba/beds/BATF3_low.bed --blacklist /home/vant/ATAC_seq/OSCC/footprinting/blacklist/hg38.blacklist.bed --outdir $bam/BATF3_low/ --cores 7

done

for bam in $(ls | egrep bam | cut -d'_' -f1,2,3,4)
do

    TOBIAS FootprintScores --signal $bam/BATF3_low/$bam*"_corrected.bw" --regions /home/vant/ATAC_seq/OSCC/prueba/beds/BATF3_low.bed \
    --output $bam/BATF3_low/$bam"_footprints_BATF3_low.bw" --cores 7

done

for bam in $(ls | egrep bam | cut -d'_' -f1,2,3,4)
do
    
    TOBIAS PlotAggregate --TFBS /home/vant/ATAC_seq/OSCC/prueba/beds/BATF3_low.bed --signals $bam/BATF3_low/$bam*"_corrected.bw" $bam/BATF3_low/$bam*"_footprints_BATF3_low.bw" --output $bam/BATF3_low/footprint.pdf --share_y both --plot_boundaries --signal-on-x --flank 250

done

fi