##Bucle ATACorrect

echo elije uno de los beds para el footprinting 

a= ls /home/scclab/Atacseq/fastq/new_htseq/Resultados/beds | cut -d '.' -f1 | cut -d '_' -f1 | uniq

echo $a 
read motif

for bam in $(ls | grep bam | cut -d'.' -f1)
do

    TOBIAS ATACorrect --bam $bam*.bam --genome /home/scclab/Atacseq/fastq/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --peaks /home/scclab/Atacseq/fastq/new_htseq/Resultados/beds/tp63_high.bed --blacklist /home/scclab/Atacseq/Bams/htseq/Resultados/motifs/blacklist/hg38-blacklist.v2.bed --outdir /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_high/ --cores 23

done

for bam in $(ls | grep bam | cut -d'.' -f1 )
do

    TOBIAS FootprintScores --signal /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_high/$bam"_corrected.bw" --regions /home/scclab/Atacseq/fastq/new_htseq/Resultados/beds/tp63_high.bed \
    --output /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_high/$bam"_footprints_tp63_high.bw" --cores 23

done

for bam in $(ls | grep bam | cut -d'.' -f1)
do
    
    TOBIAS PlotAggregate --TFBS /home/scclab/Atacseq/fastq/new_htseq/Resultados/beds/tp63_high.bed --signals /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_high/$bam*"_corrected.bw" /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_high/$bam*"_footprints_tp63_high.bw" --output /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_high/footprint.pdf --share_y both --plot_boundaries --signal-on-x --flank 250

done

for bam in $(ls | grep bam | cut -d'.' -f1)
do

    TOBIAS ATACorrect --bam $bam*.bam --genome /home/scclab/Atacseq/fastq/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --peaks /home/scclab/Atacseq/fastq/new_htseq/Resultados/beds/tp63_low.bed --blacklist /home/scclab/Atacseq/Bams/htseq/Resultados/motifs/blacklist/hg38-blacklist.v2.bed --outdir /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_low/ --cores 23

done

for bam in $(ls | grep bam | cut -d'.' -f1 )
do

    TOBIAS FootprintScores --signal /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_low/$bam"_corrected.bw" --regions /home/scclab/Atacseq/fastq/new_htseq/Resultados/beds/tp63_low.bed \
    --output /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_low/$bam"_footprints_tp63_low.bw" --cores 23

done

for bam in $(ls | grep bam | cut -d'.' -f1)
do
    
    TOBIAS PlotAggregate --TFBS /home/scclab/Atacseq/fastq/new_htseq/Resultados/beds/tp63_low.bed --signals /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_low/$bam*"_corrected.bw" /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_low/$bam*"_footprints_tp63_low.bw" --output /home/scclab/Atacseq/fastq/new_htseq/Resultados/footprinting/$bam/tp63_low/footprint.pdf --share_y both --plot_boundaries --signal-on-x --flank 250

done

