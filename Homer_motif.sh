cd /home/vant/ATAC_seq/OSCC/peak_calling/motif

for bed in $(ls | grep .bed | cut -d '.' -f1)
do
mkdir -p motif/$bed
findMotifs.pl $bed.fa human motif/$bed/ -fasta scrambled.fasta
done
fi

