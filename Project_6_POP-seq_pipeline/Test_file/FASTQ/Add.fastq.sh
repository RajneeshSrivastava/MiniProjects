
#!/bin/bash
#PBS -l nodes=1:ppn=16 
#PBS -l walltime=04:00:00 
#PBS -l gres=ccm
#PBS -N Add

module load ccmrun
module load samtools/1.9
module load gcc/5.3.0
module load hisat/1.6.2
module load bedtools/1.2
module load piranha/1.6.1

cd ./Test_file/FASTQ

ccmrun hisat -p 16 -q -x Reference/human/index/hg38* -u ./Test_file/FASTQ/Add.fastq -S Add.sam
ccmrun samtools view -bS Add.sam > Add.bam
ccmrun samtools sort Add.bam Add.sorted
ccmrun samtools index Add.sorted.bam
ccmrun bedtools bamtobed -i Add.sorted.bam > Addsorted.bed
ccmrun piranha-1.2.1/bin/Piranha -l -s Add.sorted.bed  -o Add_peaks.bed
