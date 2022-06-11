#!/bin/bash

#$ -S /bin/bash
#$ -N scov2_vars
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=1:00:0
#$ -t 1-540
#$ -j y
#$ -l tscratch=3G
#$ -o /SAN/ballouxlab/tb_lima/temp/


date
hostname

BCFTOOLS=/share/apps/genomics/bcftools-1.9/bin/bcftools;
DEV=/SAN/scov2/illumina/cluster;
REF=/SAN/scov2/reference/sarscov2.wuhan1.NC_045512.2.fa;
LOWCOMP=/SAN/scov2/reference/scov2.lowComp.bed;

IN_DIR=/SAN/scov2/illumina/bam;
OUT_DIR=/SAN/scov2/illumina/vcf;

mkdir -p $OUT_DIR;

TEMP_DIR=/scratch0/atorreso/$JOB_ID.$SGE_TASK_ID
mkdir -p $TEMP_DIR

function finish {
  rm -rf /scratch0/atorreso/$JOB_ID.$SGE_TASK_ID
}

trap finish EXIT ERR



BAM=$(ls $IN_DIR/*.bam | awk "NR==${SGE_TASK_ID}");
SAMPLE=$(basename $BAM .bam)

echo "Running mpileup"

# Call variants

$BCFTOOLS mpileup -a AD,ADF,ADR,DP,SP -A -B -x -pm3 -d 1000000 -L 1000000 -f $REF -q 30 -Q 30 $BAM | bcftools call -A -m --ploidy 1 -Oz -o $TEMP_DIR/$SAMPLE.vcf.gz

echo "Normalizing variants"
$BCFTOOLS norm -f $REF $TEMP_DIR/$SAMPLE.vcf.gz -Oz -o $TEMP_DIR/$SAMPLE.norm.vcf.gz;
tabix -p vcf $TEMP_DIR/$SAMPLE.norm.vcf.gz;

echo "Soft-filtering variants"
$BCFTOOLS filter -m+ -s'SnpGap' -g 2 $TEMP_DIR/$SAMPLE.norm.vcf.gz | bcftools filter -m+ -s'MinMQ' -e 'FMT/GT !~ "0" & INFO/MQ < 30' | bcftools filter -m+ -s'QUAL' -e 'FMT/GT !~ "0" & QUAL < 30' | bcftools filter -m+ -s'PosBias' -e 'FMT/GT !~ "0" & INFO/RPB < 0.001' | bcftools filter -m+ -s'MqBias' -e 'FMT/GT !~ "0" & INFO/MQB < 0.001' | bcftools filter -m+ -s'StrandBias' -e 'FMT/GT !~ "0" & FMT/SP > 50' | bcftools filter -m+ -s'MinDP' -e "INFO/DP < 50" -Oz -o $TEMP_DIR/$SAMPLE.sfilt.vcf.gz;

echo "Mask low complexity regions"

$DEV/softFilt_vcf_by_bed.py -i $TEMP_DIR/$SAMPLE.sfilt.vcf.gz -b $LOWCOMP -o $TEMP_DIR/$SAMPLE.sfilt.lc.vcf.gz
tabix -p vcf $TEMP_DIR/$SAMPLE.sfilt.lc.vcf.gz

echo "Add FT tag"
$DEV/addFT_pileup.py -i $TEMP_DIR/$SAMPLE.sfilt.lc.vcf.gz -o $TEMP_DIR/$SAMPLE.cons.vcf.gz -d 20 -f 0.6
$DEV/addFT_pileup.py -i $TEMP_DIR/$SAMPLE.sfilt.lc.vcf.gz -o $TEMP_DIR/$SAMPLE.hets.vcf.gz -d 20 -f 0.999


tabix -p vcf $TEMP_DIR/$SAMPLE.cons.vcf.gz
tabix -p vcf $TEMP_DIR/$SAMPLE.hets.vcf.gz


echo "Moving output files to $OUT_DIR"

mv $TEMP_DIR/$SAMPLE.cons.vcf.gz* $OUT_DIR
mv $TEMP_DIR/$SAMPLE.hets.vcf.gz* $OUT_DIR


date
hostname
