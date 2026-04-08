#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -m n
#$ -M wuw@sibcb.ac.cn
#$ -notify
#$ -N minimap2_insertion

date
PATH=$PATH:/sibcb1/wuweilab1/liangyu/Software/minimap2-2.17_x64-linux
echo $sample
dm6_FA=/sibcb1/wuweilab1/wuwei-lab/Reference/Fly/dm6/Sequence/WholeGenomeFasta/genome.fa
TE_FA=/sibcb1/wuweilab1/wuwei2022/Wuwei/WangLuLab/Drosophila_TE_full_with_HMS_GFPreporter.nochrM.fa

minimap2 -t 24 -ax map-ont $TE_FA ${sample}.pass.fq.gz > ${sample}.TE.sam
samtools view -@ 24 -F 4 -bSh  ${sample}.TE.sam > ${sample}.TE.bam
samtools fastq -@ 24 ${sample}.TE.bam > ${sample}.TE.fastq
rm ${sample}.TE.sam ${sample}.TE.bam

minimap2 -t 24 -ax map-ont $dm6_FA ${sample}.TE.fastq > ${sample}.TE.dm6.sam
samtools view -@ 24 -F 4 -bSh ${sample}.TE.dm6.sam > ${sample}.TE.dm6.bam
samtools fastq -@ 24 ${sample}.TE.dm6.bam >  ${sample}.TE.dm6.fastq
rm ${sample}.TE.dm6.sam  ${sample}.TE.dm6.bam
rm  ${sample}.TE.fastq

minimap2 -t 24 -x map-ont $TE_FA ${sample}.TE.dm6.fastq > ${sample}.TE.dm6.TE.paf
minimap2 -t 24 -x map-ont $dm6_FA ${sample}.TE.dm6.fastq > ${sample}.TE.dm6.dm6.paf

date
