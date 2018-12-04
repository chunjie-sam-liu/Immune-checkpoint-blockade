#下载cell的WES数据，SRP067938和SRP090294：
nohup prefetch --option-file 067938.ids &
nohup prefetch --option-file 090294.ids &

#WES数据转为fastq格式
cat 067938.ids |while read line;do echo "nohup fastq-dump --split-3 -gzip  --outdir /data/liucj/data/immune-checkpoint-blockade/SRP067938_fastq/  --split-files  /data/liucj/data/immune-checkpoint-blockade/SRP067938/$line"".sra &" ;done >work.sh
sh work.sh

cat 090294.ids |while read line;do echo "nohup fastq-dump --split-3 -gzip  --outdir /data/liucj/data/immune-checkpoint-blockade/SRP090294_fastq/  --split-files  /data/liucj/data/immune-checkpoint-blockade/SRP090294/$line"".sra &" ;done >work.sh
sh work.sh
