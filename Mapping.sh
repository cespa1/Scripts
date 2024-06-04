
export WD=$(pwd)
cd $WD
echo "############ Starting mapping and sorting at $(date +'%H:%M:%S')... ##############"
if [ ! -d star_index/ ]  
then 
    mkdir star_index 
fi 

if [ -f star_index/Genome ]
then    
    echo "Ya existe un index, pasando al mapping"
else
##Falta poner el codigo de la creación del index de STAR
fi

STAR --runThreadN 6 --genomeDir star_index/ --outSAMtype SAM \ 
--readFilesIn fastq/$1/*_1.fq.gz fastq/$1/*_2.fq.gz \ 
--readFilesCommand zcat --outFileNamePrefix ./out/$1/$1.sam --limitGenomeGenerateRAM 18000000000

##Removing dups with Picard 

picard -I ./out/$1/$1.sam --REMOVE_DUPLICATES=true -O ./out/$1/$1"dedup.sam" ##MIRAR EL OUTPUT

##Sorting, transformar a bam y generar el index con SAMTOOLS

samtools sort out/$1/*dedup.sam -o out/$1/$1"dedup_sorted.sam" ##MIRAR EL OUTPUT

samtools view -bh out/$1/*dedup_sorted.sam > out/$1/$1"dedup_sorted.bam" ##MIRAR EL OUTPUT

samtools index out/$1/*dedup_sorted.bam out/$1/$1"dedup_sorted.bam.bai" ##MIRAR EL OUTPUT

echo "############ Finish mapping and sorting at $(date +'%H:%M:%S')... ##############"

