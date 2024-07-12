


##Mapping y .bam .bai
cd fastq    
for muestra in $(ls | grep .fq.gz | cut -d '_' -f1,2,3,4| uniq )
do
    echo "###########################    comenzando mapping de $muestra  $(date +'%H:%M:%S')   ###########################"
    bowtie2 -x human_index -1 $muestra*trimmed_dup_1.fq.gz -2 $muestra*trimmed_dup_2.fq.gz -S $muestra".sam" --very-sensitive -X 500 -p 20 > $muestra"bowtie_log.txt"
    mv $muestra*trimmed_dup_1.fq.gz fastq_trimmed/
    mv $muestra*trimmed_dup_2.fq.gz fastq_trimmed/
    mv $muestra"bowtie_log.txt" bowtie_logs/
    echo "###########################    terminando mapping de $muestra  $(date +'%H:%M:%S')   ###########################"

    echo "###########################    comenzando creación del .bam .bai de $muestra  $(date +'%H:%M:%S')   ###########################"
    samtools sort $muestra".sam" -o $muestra"_sort.sam"   
    rm  $muestra".sam"
    picard MarkDuplicates I=$muestra"_sort.sam" O=$muestra"_sort_dedup.sam" M=picard_log.txt REMOVE_DUPLICATES=true
    rm $muestra"_sort.sam"
    samtools view -bhS $muestra"_sort_dedup.sam" > $muestra"_sort_dedup.bam"
    rm $muestra"_sort_dedup.sam"
    samtools index $muestra"_sort_dedup.bam" $muestra"_sort_dedup.bai"
    cp $muestra"_sort_dedup.bam" ../Bams/
    echo "###########################     .bam .bai hecho de $muestra  $(date +'%H:%M:%S')  ###########################"
done



##Creación del tdf

for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
    igvtools count -z 5 -w 25 -e 250 $muestra"_sort_dedup.bam" $muestra"_dedup_sorted.tdf" Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
done


##Peak_calling

for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
        macs2 callpeak -t $muestra*dedup.bam --outdir $muestra"_ATAC_peaks.broadPeak" \
    -f BAM -n ATAC -g hs -q 0.0001 --nomodel --shift -37 --extsize 73 --broad

done


cd ../Bams/
##Creación del archivo .npz para hacer el heatmap correlation

computeMatrix reference-point --referencePoint center -p 23 -S O9_1.bw O9_2.bw O31_1.bw O31_2.bw O32_1.bw O15.bw -R ../../high_DEG ../../Low_DEG  -b 3000 -a 3000 -bs=1 -p=max --outFileName test2.matrix.gz

plotHeatmap -m test2.matrix.gz --outFileName reclass2.plot.pdf --yMax 40 --colorMap=Blues

computeMatrix reference-point --referencePoint center -p 23 -S O24_1.bw O24_2.bw O35_1.bw O14_2.bw O14_3.bw O14_1.bw -R ../../high_DEG ../../Low_DEG -b 3000 -a 3000 -bs=1 -p=max --outFileName low_no_rna_test.matrix.gz

plotHeatmap -m low_no_rna_test.matrix.gz --outFileName low_no_rna.plot.pdf --yMax 40 --colorMap=Blues

computeMatrix reference-point --referencePoint center -p 23 -S O57.bw O64.bw O67.bw O55_1.bw O55_2.bw O41_1.bw O41_2.bw -R ../../high_DEG ../../Low_DEG -b 3000 -a 3000 -bs=1 -p=max --outFileName high_no_rna_test.matrix.gz

plotHeatmap -m high_no_rna_test.matrix.gz --outFileName high_no_rna.plot.pdf --yMax 40 --colorMap=Blues

computeMatrix reference-point --referencePoint center -p 23 -S O68_2.bw O25_1.bw O25_2.bw O66_1.bw O66_2.bw O82_1.bw O82_2.bw O37_1.bw -R ../../high_DEG ../../Low_DEG -b 3000 -a 3000 -bs=1 -p=max --outFileName muestras_mal_clasificadas.matrix.gz

plotHeatmap -m muestras_mal_clasificadas.matrix.gz  --outFileName muestras_mal_clasificadas.plot.pdf --yMax 40 --colorMap=Blues

computeMatrix reference-point --referencePoint center -p 20 -S ORG_11.bw ORG_18.bw ORG_37.bw ORG_43.bw -R ../../high_DEG ../../Low_DEG -b 3000 -a 3000 -bs=1 -p=max --outFileName organoides.matrix.gz

plotHeatmap -m organoides.matrix.gz --outFileName organoides.plot.pdf --yMax 40 --colorMap=Blues

computeMatrix reference-point --referencePoint center -p 20 -S O68_1.bw -R ../../high_DEG ../../Low_DEG -b 3000 -a 3000 -bs=1 -p=max --outFileName high.matrix.gzs