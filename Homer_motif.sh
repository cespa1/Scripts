cd /home/vant/ATAC_seq/OSCC/peak_calling/motif
echo "Usar scramble fasta?"
read texto

if [ $texto == no ]
then
#for bed in $(ls | grep .bed | cut -d '.' -f1)
#do
#bedtools getfasta -fi Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed $bed.bed > $bed.fa
#done


for bed in $(ls | grep .bed | cut -d '.' -f1)
do
mkdir -p motif/$bed"_no_scrambled"
findMotifs.pl $bed.fa fasta motif/$bed"_no_scrambled"/
done
else
#for bed in $(ls | grep .bed | cut -d '.' -f1)
#do
#bedtools getfasta -fi Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed $bed.bed > $bed.fa
#done


for bed in $(ls | grep .bed | cut -d '.' -f1)
do
mkdir -p motif/$bed
findMotifs.pl $bed.fa human motif/$bed/ -fasta scrambled.fasta
done
fi