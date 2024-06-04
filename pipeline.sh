## $conda activate peak_calling

##Working directory as var y creación de directorios si no existen
echo "############ Starting pipeline at $(date +'%H:%M:%S')... ##############"
cd /home/vant/ATAC_seq
export WD=$(pwd)
cd $WD
mkdir -p out/

##Preprocesamiento de los fastq, (generar achivo fastQC,trimming)
    ##DONE 
for sample in $(ls fastq/)
do
    bash /home/vant/ATAC_seq/scripts/scripts_ATAC-seq/fastq_preprocess.sh $sample
done 

##MAPPING and sorting de la muestra 
cd $WD

for sample in $(ls fastq/)
do
    bash /home/vant/ATAC_seq/scripts/scripts_ATAC-seq/scripts/Mapping.sh $sample
done

##IGVtools count
for sample in $(ls fastq/)
do
    igvtools count -z 5 -w 25 -e 250 out/$sample/$sample*edup_sorted.bam out/$sample/$sample"_dedup_sorted.tdf" human_genome/GRCh38.p13.genome.fa

done


##Peak calling con MACS2

#macs2 callpeak -t out/$muestra/*dedup.bam --outdir out/$muestra/peak_calling/ \
#   -f BAM -n ATAC -g hs -q 0.0001 --nomodel --shift -37 --extsize 73 --broad
 
##Creación de los archivos bed del peak calling y creación del peakcom

#bash /home/vant/ATAC_seq/scripts/scripts_ATAC-seq/bed_creator.sh

 echo "############ Finish pipeline at $(date +'%H:%M:%S')... ##############"