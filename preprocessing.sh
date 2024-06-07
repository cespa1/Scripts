##Entrar a la carpeta con los fastq raw
cd fastq_raw

##Control de calidad de las muestras pre-trimming
for muestra in $(ls | grep .fq.gz| cut -d '_' -f1,2,3,4| uniq)
    do
        fastqc $muestra*_1.fq.gz -o ../fastqc/ --memory 10000 --threads 20
        fastqc $muestra*_2.fq.gz -o ../fastqc/ --memory 10000 --threads 20
    done

##Eliminación de adaptadores 
for muestra in $(ls | grep .fq.gz| cut -d '_' -f1,2,3,4| uniq)
     do
        fastp -i $muestra*_1.fq.gz -o $muestra"_trimmed_dup_1.fq.gz" \
            -I $muestra*_2.fq.gz -O $muestra"_trimmed_dup_2.fq.gz" \
            --detect_adapter_for_pe --thread  16 -h $muestra".html"
        mv $muestra".html" ../fastp_reports/
        mv $muestra"_trimmed_dup_1.fq.gz" ../fastq_trimmed
        mv $muestra"_trimmed_dup_2.fq.gz" ../fastq_trimmed
    done
cd ..
cd fastq_trimmed/
##Control de calidad de las muestras post-trimming
for muestra in $(ls | grep .fq.gz| cut -d '_' -f1,2,3,4| uniq)
     do
        fastqc $muestra*"_trimmed_dup_1.fq.gz" -o ../fastqc/ --memory 10000 --threads 20
        fastqc $muestra*"_trimmed_dup_2.fq.gz" -o ../fastqc/ --memory 10000 --threads 20
    done
cd ..
##Salir al directorio de trabajo     
cd ..














