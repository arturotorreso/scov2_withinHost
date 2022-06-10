#!/bin/bash

#$ -S /bin/bash
#$ -N scov_raxml
#$ -pe mpi 6
#$ -R y
#$ -l h_vmem=2G
#$ -l tmem=2G
#$ -l h_rt=2:0:0
#$ -j y
#$ -t 1-2
#$ -o /SAN/ballouxlab/tb_lima/temp/

RAXML=/share/apps/genomics/raxml-ng-0.9.0/bin/raxml-ng-mpi;
MPIRUN=/share/apps/openmpi-3.1.1/bin/mpirun;

IN_DIR=/SAN/scov2/illumina/alignment/;
OUT_DIR=/SAN/scov2/illumina/phylogenetics/ml_tree/;

mkdir -p $OUT_DIR


# Model 1: Consensus model

if [[ $SGE_TASK_ID -eq 1 ]]; then
  echo "Consensus model"
  # $RAXML --parse --msa $IN_DIR/scov2.consensus.pass.fa --model GTR+G --prefix $IN_DIR/scov2.consensus --seed $RANDOM --precision 12 --blmin 0.000000001

  IN_FILE=$IN_DIR/scov2.consensus.raxml.rba;
  PREFIX=scov2.cons;

  $MPIRUN --map-by node -np $NSLOTS --mca btl tcp,self $RAXML --msa $IN_FILE --model GTR+G --prefix $OUT_DIR/$PREFIX --tree rand{10},pars{10} --threads 1 --seed $RANDOM --precision 12 --blmin 0.000000001;
fi;



# Model 2: Full model

if [[ $SGE_TASK_ID -eq 2 ]]; then
  echo "Full model"
  # $RAXML --parse --msa $IN_DIR/scov2.full.pass.fa --model PROTGTR+G --prefix $IN_DIR/scov2.full --seed $RANDOM --precision 12 --blmin 0.000000001

  IN_FILE=$IN_DIR/scov2.full.raxml.rba;
  PREFIX=scov2.full;

  $MPIRUN --map-by node -np $NSLOTS --mca btl tcp,self $RAXML --msa $IN_FILE --model PROTGTR+G --prefix $OUT_DIR/$PREFIX --tree pars{10},rand{10} --threads 1 --seed $RANDOM --precision 12 --blmin 0.000000001;
fi;
