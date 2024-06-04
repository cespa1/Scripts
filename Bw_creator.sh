
for muestra in $(ls | egrep bam | cut -d'_' -f1,2,3,4)
    do
        bamCoverage -b $muestra*.bam -o Bigwig/$muestra"_OSCC.bw" 
    done
