echo elije uno de los beds para el footprinting 

a= ls /home/scclab/Atacseq/Bams/htseq/Resultados/beds | cut -d '.' -f1 | cut -d '_' -f1 | uniq

echo $a 
read motif


for bam in $(ls | egrep bam | cut -d'_' -f1)
do

    TOBIAS ATACorrect --bam $bam*.bam --genome /home/scclab/Atacseq/fastq/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --peaks /home/scclab/Atacseq/Bams/htseq/Resultados/beds/$motif"_high.bed" --blacklist /home/scclab/Atacseq/Bams/htseq/Resultados/motifs/blacklist/hg38-blacklist.v2.bed --outdir /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_high"/ --cores 23

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do

    TOBIAS FootprintScores --signal /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_high"/$bam"_merged_corrected.bw" --regions /home/scclab/Atacseq/Bams/htseq/Resultados/beds/$motif"_high.bed" \
    --output /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_high"/$bam"_footprints_"$motif"_high.bw" --cores 23

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do
    
    TOBIAS PlotAggregate --TFBS /home/scclab/Atacseq/Bams/htseq/Resultados/beds/$motif"_high.bed" --signals /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_high/"$bam*"_merged_corrected.bw" /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_high"/$bam*"_footprints_"$motif"_high.bw" --output /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_high"/footprint_250.pdf --share_y both --plot_boundaries --signal-on-x --flank 250

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do

    TOBIAS ATACorrect --bam $bam*.bam --genome /home/scclab/Atacseq/fastq/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --peaks /home/scclab/Atacseq/Bams/htseq/Resultados/beds/$motif"_low.bed" --blacklist /home/scclab/Atacseq/Bams/htseq/Resultados/motifs/blacklist/hg38-blacklist.v2.bed --outdir /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_low"/ --cores 23

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do

    TOBIAS FootprintScores --signal /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_low"/$bam"_merged_corrected.bw" --regions /home/scclab/Atacseq/Bams/htseq/Resultados/beds/$motif"_low.bed" \
    --output /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_low"/$bam"_footprints_"$motif"_low.bw" --cores 23

done

for bam in $(ls | egrep bam | cut -d'_' -f1)
do
    
    TOBIAS PlotAggregate --TFBS /home/scclab/Atacseq/Bams/htseq/Resultados/beds/$motif"_low.bed" --signals /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_low/"$bam*"_merged_corrected.bw" /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_low"/$bam*"_footprints_"$motif"_low.bw" --output /home/scclab/Atacseq/Bams/htseq/Resultados/footprinting/$bam/$motif"_low"/footprint_250.pdf --share_y both --plot_boundaries --signal-on-x --flank 250

done


