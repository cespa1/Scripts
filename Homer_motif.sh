cd /home/vant/ATAC_seq/OSCC/peak_calling/motif

for bed in $(ls | grep .bed | cut -d '.' -f1)
do
mkdir -p motif/$bed
findMotifs.pl /home/scclab/Atacseq/fastq/new_htseq/Resultados/beds/low_deseq.fa human /home/scclab/Atacseq/fastq/new_htseq/Resultados/beds/homer_low/ -fasta ~/Atacseq/scrambled.fasta
done
fi

