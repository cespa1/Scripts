#1er paso, calcular footprinting de cada muestra 
    #Footprinting high samples

echo ""

cd High
for sample in $(ls | grep _high.bam | cut -d '_' -f1,2,3,4 | uniq )
do
rgt-hint footprinting --atac-seq --paired-end --organism=hg38 --output-location=./ --output-prefix=$sample"_high" $sample*".bam" $sample*".broadPeaks"
done
cd ..
    #Footprinting low samples
cd Low
for sample in $(ls | grep _low.bam | cut -d '_' -f1,2,3,4 | uniq )
do
rgt-hint footprinting --atac-seq --paired-end --organism=hg38 --output-location=./ --output-prefix=$sample"_low" $sample*".bam" $sample*".broadPeaks"
done
cd ..
echo ""

#Obtención de los footprintings por cada una de las muestras 

for sample in $(ls | grep _high.bam | cut -d '_' -f1,2,3,4 | uniq )
do
rgt-motifanalysis matching --organism=hg38 --input-files $sample*".bed"
done

#Activity score con la división en dos grupos: pEmt high y pEmt low

rgt-hint differential --organism=hg38 --bc --nc 7 --mpbs-files=./match/,--reads-files= --conditions= --output-location=High_vs_Low


