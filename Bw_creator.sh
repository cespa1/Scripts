
for muestra in $(ls | grep bam | cut -d'_' -f1,2,3,4)
    do
        bamCoverage -p 23 -b $muestra*.bam -o Bigwig/$muestra"_OSCC.bw" 
    done
