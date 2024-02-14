#!/bin/bash

#$ -S /bin/bash
#$ -N refBias_art
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=5:0:0
#$ -j y
#$ -t 1-163
#$ -l tscratch=10G
#$ -o /temp

hostname
date

TEMP_DIR=$JOB_ID.$SGE_TASK_ID;
mkdir -p $TEMP_DIR;

function finish {
  rm -rf $JOB_ID.$SGE_TASK_ID
}

trap finish EXIT ERR


REF=reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.toplevel.fa;
IN_DIR=bias_ref/assembly/ncbi-genomes-2022-05-09;
OUT_DIR=bias_ref/assembly/art_sims;

mkdir -p $OUT_DIR

SAMPLE=$(awk "NR==${SGE_TASK_ID}" bias_ref/assembly/assemblies_final.txt)
FILE=$IN_DIR/$SAMPLE/${SAMPLE}_genomic.fna.gz

zcat $FILE > $TEMP_DIR/${SAMPLE}_genomic.fna;
FILE=$TEMP_DIR/${SAMPLE}_genomic.fna;

for n in {1..10}; do
  art_illumina -ss HSXt -p -i $FILE -l 150 -f 50 -m 650 -s 150 -na -o $TEMP_DIR/$SAMPLE.rep${n};
  ls -l $TEMP_DIR;
  gzip -c $TEMP_DIR/$SAMPLE.rep${n}1.fq > $OUT_DIR/$SAMPLE.rep${n}_1.fastq.gz;
  gzip -c $TEMP_DIR/$SAMPLE.rep${n}2.fq > $OUT_DIR/$SAMPLE.rep${n}_2.fastq.gz;
done

hostname
date
