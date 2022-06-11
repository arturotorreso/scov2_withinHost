#!/bin/bash

#$ -S /bin/bash
#$ -N nano_pseudo
#$ -l h_vmem=7G
#$ -l tmem=7G
#$ -l h_rt=1:00:0
#$ -t 1-540
#$ -j y
#$ -l tscratch=3G
#$ -o /SAN/ballouxlab/tb_lima/temp/


date
hostname

REF=/SAN/scov2/reference/sarscov2.wuhan1.NC_045512.2.fa;

IN_DIR=/SAN/scov2/illumina/bam;
VCF_DIR=/SAN/scov2/illumina/vcf
OUT_MSA=/SAN/scov2/illumina/pseudosequence;


mkdir -p $OUT_MSA;

MASK_FILE=/SAN/scov2/reference/scov2_masking.txt;


TEMP_DIR=/scratch0/atorreso/$JOB_ID.$SGE_TASK_ID
mkdir -p $TEMP_DIR

function finish {
  rm -rf /scratch0/atorreso/$JOB_ID.$SGE_TASK_ID
}

trap finish EXIT ERR


BAM=$(ls $IN_DIR/*.bam | awk "NR==${SGE_TASK_ID}");
SAMPLE=$(basename $BAM .bam)


echo "Get pseudosequence for sample $SAMPLE"

/SAN/ballouxlab/tb_lima/dev/vcf2pseudoseq.py -v $VCF_DIR/$SAMPLE.cons.vcf.gz -o $TEMP_DIR/$SAMPLE.consensus.fa -H N --fail_as_N --het_freq 0.6 -d -w -r $REF -f "FORMAT/FT == 0" -f "FILTER != PASS"
/SAN/ballouxlab/tb_lima/dev/vcf2pseudoseq.py -v $VCF_DIR/$SAMPLE.hets.vcf.gz -o $TEMP_DIR/$SAMPLE.hets.fa -H iupac --fail_as_N --het_freq 0.02 -d -w -r $REF -f "FORMAT/FT == 0" -f "FILTER != PASS"


echo "Mask by depth for sample $SAMPLE"

samtools depth -aa $BAM > $TEMP_DIR/$SAMPLE.depth.txt
awk -F"\t" '$3<30' $TEMP_DIR/$SAMPLE.depth.txt | cut -f2 > $TEMP_DIR/pos.unsrt.txt

sort -nu $TEMP_DIR/pos.unsrt.txt > $TEMP_DIR/pos.txt

seqkit mutate $(while read pos; do echo "-p ${pos}:N"; done <<< "$(sort -n $TEMP_DIR/pos.txt)" | tr '\n' ' ') -w 0 $TEMP_DIR/$SAMPLE.consensus.fa --quiet -o $TEMP_DIR/$SAMPLE.consensus.masked.fa;
seqkit mutate $(while read pos; do echo "-p ${pos}:N"; done <<< "$(sort -n $TEMP_DIR/pos.txt)" | tr '\n' ' ') -w 0 $TEMP_DIR/$SAMPLE.hets.fa --quiet -o $TEMP_DIR/$SAMPLE.hets.masked.fa;


echo "Adding major/minor variant information"

bcftools view -i'FMT/GT !~ "0/0"' $VCF_DIR/$SAMPLE.hets.vcf.gz -Oz -o $TEMP_DIR/$SAMPLE.vars.vcf.gz
/SAN/ballouxlab/tb_lima/dev/add_maj_min.py -i $TEMP_DIR/$SAMPLE.hets.masked.fa -v $TEMP_DIR/$SAMPLE.vars.vcf.gz -o $TEMP_DIR/$SAMPLE.full.masked.fa;


echo "Moving output files for sample $SAMPLE"

mv $TEMP_DIR/$SAMPLE.consensus.masked.fa $OUT_MSA/$SAMPLE.consensus.fa;
mv $TEMP_DIR/$SAMPLE.hets.masked.fa $OUT_MSA/$SAMPLE.hets.fa
mv $TEMP_DIR/$SAMPLE.full.masked.fa $OUT_MSA/$SAMPLE.full.fa

date
hostname
