##Bucle ATACorrect

for bam in $(ls | egrep bam | cut -d'_' -f1)
do

    TOBIAS ATACorrect --bam $bam*.bam --genome /home/scclab/Atacseq/fastq/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --peaks /home/scclab/Atacseq/Bams/htseq/Resultados/beds/PRDM1_high.bed --blacklist /home/scclab/Atacseq/Bams/htseq/Resultados/motifs/blacklist/hg38-blacklist.v2.bed --outdir /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_high/ --cores 23

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do

    TOBIAS FootprintScores --signal /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_high/$bam"_merged_corrected.bw" --regions /home/scclab/Atacseq/Bams/htseq/Resultados/beds/PRDM1_high.bed \
    --output /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_high/$bam"_footprints_PRDM1_high.bw" --cores 23

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do
    
    TOBIAS PlotAggregate --TFBS /home/scclab/Atacseq/Bams/htseq/Resultados/beds/PRDM1_high.bed --signals /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_high/$bam*"_merged_corrected.bw" /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_high/$bam*"_footprints_PRDM1_high.bw" --output /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_high/footprint.pdf --share_y both --plot_boundaries --signal-on-x --flank 250

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do

    TOBIAS ATACorrect --bam $bam*.bam --genome /home/scclab/Atacseq/fastq/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --peaks /home/scclab/Atacseq/Bams/htseq/Resultados/beds/PRDM1_low.bed --blacklist /home/scclab/Atacseq/Bams/htseq/Resultados/motifs/blacklist/hg38-blacklist.v2.bed --outdir /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_low/ --cores 23

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do

    TOBIAS FootprintScores --signal /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_low/$bam*"_merged_corrected.bw" --regions /home/scclab/Atacseq/Bams/htseq/Resultados/beds/PRDM1_low.bed \
    --output /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_low/$bam"_footprints_PRDM1_low.bw" --cores 2

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do
    
    TOBIAS PlotAggregate --TFBS /home/scclab/Atacseq/Bams/htseq/Resultados/beds/PRDM1_low.bed --signals /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_low/$bam*"_merged_corrected.bw" /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_low/$bam*"_footprints_PRDM1_low.bw" --output /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/PRDM1_low/footprint.pdf --share_y both --plot_boundaries --signal-on-x --flank 250

done