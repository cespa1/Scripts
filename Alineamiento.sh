
##Creación del index si no existe previamente 
cd index_DB
if [ -f fastq/human_index.1.bt2 ]

then 
   echo "Index humano ya existe"

else
echo "###########################    Haciendo index de humano   ###########################"
curl -O http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 
bowtie2-build Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz human_index
fi
cd ..

##Alineamiento de muestras y creación de los .bam y .bai

cd fastq_trimmed
for muestra in $(ls | grep .fq.gz | cut -d '_' -f1,2,3,4| uniq )
do
    echo "###########################    comenzando mapping de $muestra  $(date +'%H:%M:%S')   ###########################"
    bowtie2 -x human_index -1 $muestra*trimmed_dup_1.fq.gz -2 $muestra*trimmed_dup_2.fq.gz -S $muestra".sam" --very-sensitive -X 500 -p 20 > $muestra"bowtie_log.txt"
    mv $muestra"bowtie_log.txt" ../bowtie_logs/
    echo "###########################    terminando mapping de $muestra  $(date +'%H:%M:%S')   ###########################"

    echo "###########################    comenzando creación del .bam .bai de $muestra  $(date +'%H:%M:%S')   ###########################"
    samtools sort $muestra".sam" -o $muestra"_sort.sam"   
    rm  $muestra".sam"
    picard MarkDuplicates I=$muestra"_sort.sam" O=$muestra"_sort_dedup.sam" M=picard_log.txt REMOVE_DUPLICATES=true
    rm $muestra"_sort.sam"
    samtools view -bh $muestra"_sort_dedup.sam" > $muestra"_sort_dedup.bam"
    rm $muestra"_sort_dedup.sam"
    samtools index $muestra"_sort_dedup.bam" $muestra"_sort_dedup.bai"
    mv $muestra"_sort_dedup.bam" ../Bam/
    mv $muestra"_sort_dedup.bai" ../Bam/
    echo "###########################     .bam .bai hecho de $muestra  $(date +'%H:%M:%S')  ###########################"
done
cd ..