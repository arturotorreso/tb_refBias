#!/bin/bash

#$ -S /bin/bash
#$ -N refBias_mapping
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=5:0:0
#$ -j y
#$ -t 1-20
#$ -l tscratch=10G
#$ -o temp

hostname
date

REF=reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.toplevel.fa;

IN_DIR=bias_ref/raw_reads/fastq;
OUT_DIR=bias_ref/raw_reads/1_mapping/vcf;

mkdir -p $OUT_DIR

ADAPTERS=software/trimmomatic/trimmomatic.fa;
MASK_FILE=reference/mtb.h37rv.hypervariable.bed;
DR_FILE=reference/tb_profiler_db.complete.bed;

GATK=software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
BEDTOOLS=software/bedtools-2.25.0/bin/bedtools;
PICARD=software/picard-2.21.9/picard.jar;

TEMP_DIR=$JOB_ID.$SGE_TASK_ID;
mkdir -p $TEMP_DIR;

function finish {
  rm -rf $JOB_ID.$SGE_TASK_ID
}

trap finish EXIT ERR


###################
#    GET FASTQ    #
###################

SAMPLE=$(basename $(ls $IN_DIR/*_1.fastq.gz | awk "NR==${SGE_TASK_ID}") _1.fastq.gz)

FWD=$(ls $IN_DIR/${SAMPLE}_1.fastq.gz)
REV=$(ls $IN_DIR/${SAMPLE}_2.fastq.gz)

trimmomatic PE $FWD $REV $TEMP_DIR/${SAMPLE}_R1.fastq.gz /dev/null $TEMP_DIR/${SAMPLE}_R2.fastq.gz /dev/null -phred33 ILLUMINACLIP:$ADAPTERS:1:30:11 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 AVGQUAL:20

FWD=$TEMP_DIR/${SAMPLE}_R1.fastq.gz;
REV=$TEMP_DIR/${SAMPLE}_R2.fastq.gz;


##################################
# Mapp reads to sample-reference #
##################################

echo "Mapping reads into sample-specific reference"

HEADER=$(zcat $FWD | head -n1 | cut -d' ' -f2-);

FLOWCELL=$(echo $HEADER | cut -d'-' -f1);

rg_id=$FLOWCELL;
pu=$FLOWCELL.$SAMPLE;
pl=ART_ILLUMINA;
sm=$SAMPLE;


echo "Aligning raw reads back to $REF"
bwa mem -M -t 2 -v 1 -R "@RG\tID:$rg_id\tPU:$pu\tPL:$pl\tSM:$sm" $REF $FWD $REV | samtools view -bS - | samtools sort -@ 4 -T $TEMP_DIR/$SAMPLE.temp -O bam -o $TEMP_DIR/$SAMPLE.raw.bam;

echo "Marking duplicates"
sambamba markdup --tmpdir=$TEMP_DIR $TEMP_DIR/$SAMPLE.raw.bam $TEMP_DIR/$SAMPLE.markDup.bam --overflow-list-size 600000
samtools index $TEMP_DIR/$SAMPLE.markDup.bam;

echo "Indel re-alignment"
java -jar $GATK -T RealignerTargetCreator -nt 4 -R $REF -I $TEMP_DIR/$SAMPLE.markDup.bam -o $TEMP_DIR/$SAMPLE.IndelRealigner.intervals;
java -jar $GATK -T IndelRealigner -R $REF -I $TEMP_DIR/$SAMPLE.markDup.bam -targetIntervals $TEMP_DIR/$SAMPLE.IndelRealigner.intervals -o $TEMP_DIR/$SAMPLE.markDup.reAlign.bam;

samtools index $TEMP_DIR/$SAMPLE.markDup.reAlign.bam;


####################
# Calling variants #
####################

BAM=$TEMP_DIR/$SAMPLE.markDup.reAlign.bam;

echo "Calling variants"

bcftools mpileup -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP -A -x -E -pm5 -d 10000 -L 10000 -f $REF -q 0 -Q 0 -Ou $BAM | bcftools call -mv -M -A -f GQ,GP -Oz | bcftools norm -f $REF -Oz | bcftools filter -m+ -s'MinMQ' -e 'FMT/GT !~ "0/0" & INFO/MQ < 20' | bcftools filter -m+ -s'QUAL' -e 'FMT/GT !~ "0/0" & QUAL < 20' | bcftools filter -m+ -s'minAD' -e 'FMT/GT ~ "1/1" & FMT/AD[:1] < 10' | bcftools filter -m+ -s'minADF' -e 'FMT/GT ~ "1/1" & FMT/ADF[:1] < 3' | bcftools filter -m+ -s'minADR' -e 'FMT/GT ~ "1/1" & FMT/ADR[:1] < 3' | bcftools filter -m+ -s'minIDV' -e 'INFO/IDV < 10' | bcftools filter -m+ -s'MinDP' -e "INFO/DP < 20" -Oz -o $TEMP_DIR/$SAMPLE.sfilt.vcf.gz;


echo "Moving output to final folder"
scp $TEMP_DIR/$SAMPLE.sfilt.vcf.gz $OUT_DIR/$SAMPLE.vcf.gz

date
