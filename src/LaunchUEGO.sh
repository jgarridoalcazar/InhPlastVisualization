#!/bin/sh
#$ -S /bin/sh
#$ -N uego-2-param
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/src/UEGOv1.0/results
#$ -j y
#$ -cwd
#$ -V
#$ -v OMP_NUM_THREADS=1
#$ -q NOParalela

./uego optiParamCereb