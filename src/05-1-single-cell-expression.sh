salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i Homo_sapiens_index

for i in $(seq 56 67);
do 

salmon quant -i Homo_sapiens_index -l A \
-r /data/liucj/data/immune-checkpoint-blockade/fastq/ERP105867_fastq/RNA-Seq/ERR22287${i}.fastq.gz \
-p 8 -o Homo_quants/ERR22287${i}_quant
done

