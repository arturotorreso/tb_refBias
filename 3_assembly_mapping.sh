#!/bin/bash

#$ -S /bin/bash
#$ -N refBias_ass_map
#$ -l h_vmem=15G
#$ -l tmem=15G
#$ -l h_rt=5:0:0
#$ -j y
#$ -t 1-1630
#$ -l tscratch=10G
#$ -o temp

hostname
date

REF=reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.toplevel.fa;

IN_DIR=bias_ref/assembly/art_sims;
OUT_DIR=bias_ref/assembly/2_assembly_mapping/vcf;
CONTIG_DIR=bias_ref/assembly/2_assembly_mapping/contigs;
BAM_DIR=bias_ref/assembly/2_assembly_mapping/bam;


mkdir -p $OUT_DIR

ADAPTERS=software/trimmomatic/trimmomatic.fa;
MASK_FILE=reference/mtb.h37rv.hypervariable.bed;
DR_FILE=reference/tb_profiler_db.complete.bed;

GATK=software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
BEDTOOLS=software/bedtools-2.25.0/bin/bedtools;


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


##################
# First assembly #
##################


echo "Assembly for $SAMPLE using SPADES"

KMER=21,33,55,77,101

echo "Running spades for $SAMPLE"
spades.py -o $TEMP_DIR -1 $FWD -2 $REV --isolate -t 2 -m 250 -k $KMER;

echo "Moving contigs from $TEMP_DIR to $CONTIG_DIR"
scp $TEMP_DIR/contigs.fasta $CONTIG_DIR/$SAMPLE.contigs.fasta;


#####################
# Filtering contigs #
#####################

echo "Filtering contigs for $SAMPLE"
CONTIGS=$TEMP_DIR/contigs.fasta

seqkit fx2tab $CONTIGS | csvtk mutate -H -t -f 1 -p "cov_(.+)" | csvtk mutate -H -t -f 1 -p "length_([0-9]+)" | awk -F "\t" '$4>=2 && $5>=300' | seqkit tab2fx > $TEMP_DIR/$SAMPLE.filt_contigs.fa

scp $TEMP_DIR/$SAMPLE.filt_contigs.fa $CONTIG_DIR/$SAMPLE.filt_contigs.fa

CONTIGS=$TEMP_DIR/$SAMPLE.filt_contigs.fa


###################
# Mapping contigs #
###################

echo "Mapping contigs to H37Rv for $SAMPLE"

minimap2 -ax asm20 -R "@RG\tID:$SAMPLE\tPU:$SAMPLE\tPL:ILLUMINA\tCN:Sanger\tSM:$SAMPLE" $REF $CONTIGS | samtools view -bS - | samtools sort -@ 4 -T $TEMP_DIR/$SAMPLE.temp -O bam -o $TEMP_DIR/$SAMPLE.algn_contig.bam;


#############################
# Sample-specific reference #
#############################

echo "Creating $SAMPLE specific reference with pseudosequence"

bcftools mpileup -a AD,ADF,ADR,DP,SP -A -x -E -pm3 -d 1000 -L 1000 -f $REF -q 0 -Q 0 -Ou $TEMP_DIR/$SAMPLE.algn_contig.bam | bcftools call -mv -M -A -f GQ,GP -Ov | vcf2pseudoseq.py -v - -o $TEMP_DIR/$SAMPLE.pseudo.contig.fa -H max --fail_as_ref -w -r $REF;

# Get Chromosome into header
sed -i "s/>/&Chromosome /" $TEMP_DIR/$SAMPLE.pseudo.contig.fa;

# Just in case, change any non-IUPAC into N
cat <(head -n1 $TEMP_DIR/$SAMPLE.pseudo.contig.fa) <(tail -n+2 $TEMP_DIR/$SAMPLE.pseudo.contig.fa | tr 'MRWSYKVHDB-' 'N') > $TEMP_DIR/$SAMPLE.pseudo.contig.fa.temp;
mv $TEMP_DIR/$SAMPLE.pseudo.contig.fa.temp $TEMP_DIR/$SAMPLE.pseudo.contig.fa;


# ##################################
# # Mapp reads to sample-reference #
# ##################################

echo "Mapping reads into sample-specific reference"

bwa index $TEMP_DIR/$SAMPLE.pseudo.contig.fa;
java -jar $PICARD CreateSequenceDictionary R=$TEMP_DIR/$SAMPLE.pseudo.contig.fa O=$TEMP_DIR/$SAMPLE.pseudo.contig.dict
samtools faidx $TEMP_DIR/$SAMPLE.pseudo.contig.fa;

HEADER=$(zcat $FWD | head -n1 | cut -d' ' -f2-);

FLOWCELL=$(echo $HEADER | cut -d'-' -f1);

rg_id=$FLOWCELL;
pu=$FLOWCELL.$SAMPLE;
pl=ART_ILLUMINA;
sm=$SAMPLE;

echo "Aligning raw reads back to $SAMPLE specific reference"
bwa mem -M -t 2 -v 1 -R "@RG\tID:$rg_id\tPU:$pu\tPL:$pl\tSM:$sm" $TEMP_DIR/$SAMPLE.pseudo.contig.fa $FWD $REV | samtools view -bS - | samtools sort -@ 4 -T $TEMP_DIR/$SAMPLE.temp -O bam -o $TEMP_DIR/$SAMPLE.toContig.bam;

echo "Marking duplicates"

sambamba markdup --tmpdir=$TEMP_DIR $TEMP_DIR/$SAMPLE.toContig.bam $TEMP_DIR/$SAMPLE.toContig.markDup.bam --overflow-list-size 600000
samtools index $TEMP_DIR/$SAMPLE.toContig.markDup.bam;

echo "Indel re-alignment"

java -jar $GATK -T RealignerTargetCreator -nt 4 -R $TEMP_DIR/$SAMPLE.pseudo.contig.fa -I $TEMP_DIR/$SAMPLE.toContig.markDup.bam -o $TEMP_DIR/$SAMPLE.IndelRealigner.intervals;
java -jar $GATK -T IndelRealigner -R $TEMP_DIR/$SAMPLE.pseudo.contig.fa -I $TEMP_DIR/$SAMPLE.toContig.markDup.bam -targetIntervals $TEMP_DIR/$SAMPLE.IndelRealigner.intervals -o $TEMP_DIR/$SAMPLE.toContig.markDup.reAlign.bam;

samtools index $TEMP_DIR/$SAMPLE.toContig.markDup.reAlign.bam;

scp $TEMP_DIR/$SAMPLE.toContig.markDup.reAlign.bam /SAN/ballouxlab/tb_lima/bias_ref/assembly/2_assembly_mapping/bam/$SAMPLE.toContig.bam
scp $TEMP_DIR/$SAMPLE.toContig.markDup.reAlign.bam.bai /SAN/ballouxlab/tb_lima/bias_ref/assembly/2_assembly_mapping/bam/$SAMPLE.toContig.bam.bai


# ####################
# # Calling variants #
# ####################

BAM=$TEMP_DIR/$SAMPLE.toContig.markDup.bam;

echo "Calling variants"

bcftools mpileup -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP -A -x -E -pm5 -d 10000 -L 10000 -f $TEMP_DIR/$SAMPLE.pseudo.contig.fa -q 0 -Q 0 -Ou $BAM | bcftools call -mv -M -A -f GQ,GP -Oz | bcftools norm -f $TEMP_DIR/$SAMPLE.pseudo.contig.fa -Oz | bcftools filter -m+ -s'MinMQ' -e 'FMT/GT !~ "0/0" & INFO/MQ < 20' | bcftools filter -m+ -s'QUAL' -e 'FMT/GT !~ "0/0" & QUAL < 20' | bcftools filter -m+ -s'minAD' -e 'FMT/GT ~ "1/1" & FMT/AD[:1] < 10' | bcftools filter -m+ -s'minADF' -e 'FMT/GT ~ "1/1" & FMT/ADF[:1] < 3' | bcftools filter -m+ -s'minADR' -e 'FMT/GT ~ "1/1" & FMT/ADR[:1] < 3' | bcftools filter -m+ -s'minIDV' -e 'INFO/IDV < 10' | bcftools filter -m+ -s'MinDP' -e "INFO/DP < 20" -Oz -o $TEMP_DIR/$SAMPLE.sfilt.vcf.gz;


echo "Moving output to final folder"
scp $TEMP_DIR/$SAMPLE.sfilt.vcf.gz $OUT_DIR/$SAMPLE.vcf.gz

date
