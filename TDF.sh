
##Creación de los archivos de visualización de IGV, TDF
cd Bam/
for muestra in $(ls | grep _sort_dedup.bam | cut -d '_' -f1,2,3,4)
do
    igvtools count -z 5 -w 25 -e 250 $muestra"_sort_dedup.bam" $muestra"_dedup_sorted.tdf" ../Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
    mv $muestra"_dedup_sorted.tdf" TDF/
done
cd ..