#!/bin/bash

#$ -S /bin/bash
#$ -N covid_align
#$ -l h_vmem=7G
#$ -l tmem=7G
#$ -l h_rt=200:0:0
#$ -j y
#$ -t 2-541
#$ -l tscratch=5G
#$ -o /SAN/ballouxlab/tb_lima/temp

date
hostname


IN_DIR=/SAN/scov2/illumina/fastq;
OUT_DIR=/SAN/scov2/illumina/bam;
OUT_QC=/SAN/scov2/illumina/qc_fastqc/post_qc;
OUT_STATS=/SAN/scov2/illumina/aln_stats;

REF=/SAN/scov2/reference/sarscov2.wuhan1.NC_045512.2.withHuman.fa;

PRIMERS=/SAN/scov2/reference/nCov2019.primerV3.artic.bed;
META=/SAN/scov2/illumina/metadata/sequence_info.csv

mkdir -p $OUT_DIR;

TEMP_DIR=/scratch0/atorreso/$JOB_ID.$SGE_TASK_ID;
mkdir -p $TEMP_DIR;

function finish {
  rm -rf /scratch0/atorreso/$JOB_ID.$SGE_TASK_ID
}

trap finish EXIT ERR


SAMPLE=$(cut -d, -f1 $META | awk "NR==${SGE_TASK_ID}");
SEQ_ID=$(cut -d, -f2 $META | awk "NR==${SGE_TASK_ID}");
LIBRARY=$(cut -d, -f5 $META | awk "NR==${SGE_TASK_ID}");

# Get all the files for these run
find $IN_DIR -name '*.fastq.gz' | grep "$SEQ_ID" | grep "R1" > $TEMP_DIR/reads.$JOB_ID.$SGE_TASK_ID.txt

# Align all the fastq files separately
while read line; do
  fastq_file=$(basename $line _R1_001.fastq.gz);
  DIR=$(dirname $line)
  fastq1=$DIR/${fastq_file}_R1_001.fastq.gz;
  fastq2=$DIR/${fastq_file}_R2_001.fastq.gz;

  FLOWCELL=$(zcat $fastq1 | head -n1 | cut -d':' -f3);
  LANE=$(zcat $fastq1 | head -n1 | cut -d':' -f4);

  rg_id=$FLOWCELL.$LANE;
  lb=$LIBRARY;
  pu=$FLOWCELL.$LANE.$SAMPLE.$SEQ_ID;
  pl=ILLUMINA;
  cn=PGU

  # Quality trim of the reads
  trimmomatic PE $fastq1 $fastq2 $TEMP_DIR/$fastq_file.forward_paired.fq.gz $TEMP_DIR/$fastq_file.forward_unpaired.fq.gz $TEMP_DIR/$fastq_file.reverse_paired.fq.gz $TEMP_DIR/$fastq_file.reverse_unpaired.fq.gz SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 AVGQUAL:20

  fastq1=$TEMP_DIR/$fastq_file.forward_paired.fq.gz;
  fastq2=$TEMP_DIR/$fastq_file.reverse_paired.fq.gz;

  bwa mem -M -t 2 -v 1 -R "@RG\tID:$rg_id\tPU:$pu\tPL:$pl\tCN:$cn\tLB:$lb\tSM:$SAMPLE.$SEQ_ID" $REF $fastq1 $fastq2 | samtools view -bS - | samtools sort -@ 4 -T $TEMP_DIR/$fastq_file.temp -O bam -o $TEMP_DIR/$fastq_file.wgs.raw.bam;
  samtools index $TEMP_DIR/$fastq_file.wgs.raw.bam;

  # Get only the virus reads (NC_045512.2)
  # -f 2 to get only reads in proper pair
  # -F 256 to exclude those that are not the primary alignment
  samtools view -b -f 2 -F 256 $TEMP_DIR/$fastq_file.wgs.raw.bam -o $TEMP_DIR/$fastq_file.wgs.raw.noHost.bam NC_045512.2;
  samtools index $TEMP_DIR/$fastq_file.wgs.raw.noHost.bam;

  # Get mapping stats
  human_rds=$(samtools view $TEMP_DIR/$fastq_file.wgs.raw.bam | cut -f3 | sort | egrep -v "NC_045512.2|\*" | wc -l);
  scov2_rds=$(samtools view $TEMP_DIR/$fastq_file.wgs.raw.bam | cut -f3 | sort | egrep "NC_045512.2" | wc -l);
  unmapped=$(samtools view $TEMP_DIR/$fastq_file.wgs.raw.bam | cut -f3 | sort | egrep "\*" | wc -l);

  printf "sample\tunmapped\tscov2\thuman\n"$SAMPLE.$SEQ_ID"\t"$unmapped"\t"$scov2_rds"\t"$human_rds"\n" > $OUT_STATS/$SAMPLE.$SEQ_ID.mapStats.txt

  # Re-header to remove all the human contigs
  cat <(samtools view -H $TEMP_DIR/$fastq_file.wgs.raw.noHost.bam | grep -Ev '@SQ.*SN:chr[0-9]|Un|chrX|chrY|chrM|HLA|chrEBV|NT_|NW_|NC_000|NC_012') <(samtools view $TEMP_DIR/$fastq_file.wgs.raw.noHost.bam) | samtools view -b - -o $TEMP_DIR/$fastq_file.wgs.raw.bam;
  samtools index $TEMP_DIR/$fastq_file.wgs.raw.bam

done < $TEMP_DIR/reads.$JOB_ID.$SGE_TASK_ID.txt


samtools merge - $TEMP_DIR/*.wgs.raw.bam | samtools view -bS - | samtools sort -@ 4 -T $TEMP_DIR/$SAMPLE.$SEQ_ID.temp -O bam -o $TEMP_DIR/$SAMPLE.$SEQ_ID.bam;

BAM=$TEMP_DIR/$SAMPLE.$SEQ_ID.bam;
samtools index $BAM;


ivar trim -b $PRIMERS -p $TEMP_DIR/$SAMPLE.$SEQ_ID.unsrt -i $BAM -q 15 -m 20 -s 4 -e
samtools view -bS $TEMP_DIR/$SAMPLE.$SEQ_ID.unsrt.bam | samtools sort -T $TEMP_DIR -O bam -o $TEMP_DIR/$SAMPLE.$SEQ_ID.bam
samtools index $TEMP_DIR/$SAMPLE.$SEQ_ID.bam;


echo "Moving files from $SAMPLE.$SEQ_ID.bam and $SAMPLE.$SEQ_ID.bam.bai from $TEMP_DIR to $OUT_DIR"

mv $TEMP_DIR/$SAMPLE.$SEQ_ID.bam* $OUT_DIR


date
hostname
