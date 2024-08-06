##Low peaks

for Bams in $(ls | grep merged.bam | cut -d '_' -f1)
do

mkdir -p nucleosomas/$Bams"_low" 
nucleoatac run --bed /home/scclab/Atacseq/fastq/new_htseq/Resultados/low_deseq2 --bam /home/scclab/Atacseq/Bams/merge_bams/$Bams"_merged.bam" --fasta /home/scclab/Atacseq/fastq/Homo_sapiens.GRCh38.dna.primary_assembly.fa --cores 20 --out ./nucleosomas/$Bams"_low"/$Bams"_low"

done

##High peaks

for Bams in $(ls | grep merged.bam | cut -d '_' -f1)
do

mkdir -p nucleosomas/$Bams"_high" 
nucleoatac run --bed /home/scclab/Atacseq/fastq/new_htseq/Resultados/high_deseq2.bed --bam /home/scclab/Atacseq/Bams/merge_bams/$Bams"_merged.bam" --fasta /home/scclab/Atacseq/fastq/Homo_sapiens.GRCh38.dna.primary_assembly.fa --cores 20 --out ./nucleosomas/$Bams"_high"/$Bams"_high" 

done



