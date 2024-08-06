for muestra in $( ls | grep .bam  |cut -d '_' -f1 | uniq )
    do
        #TOBIAS ATACorrect --bam $muestra"_merged.bam" --genome  --peaks /home/scclab/Atacseq/fastq/new_htseq/Resultados/high_deseq2 --blacklist /home/scclab/Atacseq/Bams/htseq/Resultados/motifs/blacklist/hg38-blacklist.v2.bed --outdir ./all_peaks/$muestra/$muestra --cores 24

        TOBIAS FootprintScores --signal ./all_peaks/$muestra/$muestra/$muestra"_merged_corrected.bw" --regions /home/scclab/Atacseq/fastq/new_htseq/Resultados/high_deseq2 \
        --output ./all_peaks/$muestra/$muestra/$muestra"_footprints.bw" --cores 23
    done

    TOBIAS BINDetect --motifs /home/scclab/Descargas/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt --signals /home/scclab/Atacseq/Bams/merge_bams/all_peaks/43T/43T_footprints.bw /home/scclab/Atacseq/Bams/merge_bams/all_peaks/O14/O14_footprints.bw \
    --genome /home/scclab/Atacseq/fastq/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --peaks /home/scclab/Atacseq/fastq/new_htseq/Resultados/high_deseq2 --outdir /home/scclab/Atacseq/Bams/merge_bams/all_peaks/BINDetect_output_43_vs_14/ --cond_names high low --cores 23