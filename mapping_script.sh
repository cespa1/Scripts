echo "¿Qué tipo de muestra estas usando raton o humano?"
read nombre_muestra
until [ $nombre_muestra == "humano" ] || [ $nombre_muestra ==  "raton" ]
do 
    echo "por favor introduce "raton" o "humano" según tu tipo de muestra"
    read nombre_muestra
done

cd fastq
echo "###########################    Empezando el preprocesamiento de los fastq    ###########################"
   for muestra in $(ls | grep .fq.gz| cut -d '_' -f1,2,3,4| uniq)
       do
           fastqc $muestra*_1.fq.gz -o fastqc/ --memory 10000 --threads 20
	       fastqc $muestra*_2.fq.gz -o fastqc/ --memory 10000 --threads 20
           fastp -i $muestra*_1.fq.gz -o $muestra"_trimmed_dup_1.fq.gz" \
            -I $muestra*_2.fq.gz -O $muestra"_trimmed_dup_2.fq.gz" \
            --detect_adapter_for_pe --thread  16 -h $muestra".html"
            mv $muestra".html" fastp_reports/
            mv $muestra*L*_*.fq.gz fq_raw/
       done
cd ..
echo "###########################    Preprocesamiento hecho    ###########################"
##Index create
if [ $nombre_muestra == "humano" ]
then
if [ -f fastq/human_index.1.bt2 ]

then 
   echo "Index humano ya existe"

else
echo "###########################    Haciendo index de humano   ###########################"
cd fastq
curl -O http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 
bowtie2-build Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz human_index
cd ..
fi
elif [ -f fastq/mouse_index.1.bt2 ]
then
echo "Index ratón ya existe"
else
echo "###########################    Haciendo index de ratón    ###########################"
cd fastq
curl -O http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
bowtie2-build MGRCm38.primary_assembly.genome.fa.gz mouse_index
cd ..
fi

cd fastq


##Mapping y .bam .bai
if [ $nombre_muestra == "humano" ]
then
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
    samtools view -bh $muestra"_sort_dedup.sam" > $muestra"_sort_dedup.bam"
    rm $muestra"_sort_dedup.sam"
    samtools index $muestra"_sort_dedup.bam" $muestra"_sort_dedup.bai"
    echo "###########################     .bam .bai hecho de $muestra  $(date +'%H:%M:%S')  ###########################"
done

else  
for muestra in $(ls | grep .fq.gz | cut -d '_' -f1,2,3,4| uniq )
do
    echo "###########################    comenzando mapping de $muestra  $(date +'%H:%M:%S')   ###########################"
    bowtie2 -x mouse_index -1 $muestra*trimmed_dup_1.fq.gz -2 $muestra*trimmed_dup_2.fq.gz -S $muestra".sam" --very-sensitive -X 2000 -p 20
    echo "###########################    terminando mapping de $muestra  $(date +'%H:%M:%S')   ###########################"

    echo "###########################    comenzando creación del .bam .bai de $muestra  $(date +'%H:%M:%S')  ###########################"
    samtools sort $muestra".sam" -o $muestra"_sort.sam"   
    rm  $muestra".sam"
    picard MarkDuplicates I=$muestra"_sort.sam" O=$muestra"_sort_dedup.sam" M=picard_log.txt REMOVE_DUPLICATES=true
    rm $archivo"_sort.sam"
    samtools view -bh $muestra"_sort_dedup.sam" > $muestra"_sort_dedup.bam"
    rm $archivo"_sort_dedup.sam"
    samtools index $muestra"_sort_dedup.bam" $muestra"_sort_dedup.bai"
    echo "###########################     .bam .bai hecho de $muestra  $(date +'%H:%M:%S') ###########################"
done   
fi

##Creación del tdf
if [ $nombre_muestra == "humano" ]
then
for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
    igvtools count -z 5 -w 25 -e 250 $muestra"_sort_dedup.bam" $muestra"_dedup_sorted.tdf" Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
done
else
for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
    igvtools count -z 5 -w 25 -e 250 $muestra*dedup_sorted.bam $muestra"_dedup_sorted.tdf" Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
done
fi


##Peak_calling
if [ $nombre_muestra == "humano" ]
then
for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
        macs2 callpeak -t $muestra*dedup.bam --outdir $muestra"_ATAC_peaks.broadPeak" \
    -f BAM -n ATAC -g hs -q 0.0001 --nomodel --bw 1000 --shift -37 --extsize 73 --broad

        macs2 callpeak -t $muestra*dedup.bam --outdir $muestra"_ATAC_peaks.broadPeak" \
    -f BAM -n ATAC -g hs -q 0.0001 --nomodel --bw 1000 --shift -37 --extsize 73 
done
else
for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
    macs2 callpeak -t $muestra*dedup.bam --outdir $muestra"_ATAC_peaks.broadPeak" \
    -f BAM -n ATAC -g mm -q 0.0001 --nomodel --shift -37 --extsize 73 --broad
done
fi

##Peakome

touch merge.bed
for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
cut -f1-4 $muestra"_ATAC_peaks.broadPeak"/ATAC_peaks.broadPeak > $muestra"_ATAC_peaks.broadPeak"/$muestra"_peaks.bed"
cat $muestra"_ATAC_peaks.broadPeak"/$muestra"_peaks.bed" >> merge.bed
done

sort -k1,1 -k2,2n merge.bed > sorted_merge.bed

bedtools merge -i sorted_merge.bed > merge_final.bed
rm merge.bed

awk '{print $1"\tunknown\texon\t"$2"\t"$3"\t.\t+\t0\tgene_id \42"$1":"$2"-"$3"\42;\ttranscript_id \42"$1":"$2"-"$3"\42;\tgene_name \42"$1":"$2"-"$3"\42;"}' merge_final.bed  > all_peaks_merged.gtf

for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do 
 htseq-count -n 20 --format=bam --type=exon --idattr=gene_name $muestra*.bam all_peaks_merged.gtf > $muestra"__htseq_counts"
done


