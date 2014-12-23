#!/bin/sh
#$ -S /bin/sh
#$ -N uego-2-param
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/src/UEGOv1.0/results
#$ -j y
#$ -cwd
#$ -V
#$ -v OMP_NUM_THREADS=32
#$ -q 120Hbigmem
#$ -pe openmpi 32

./uego optiParamCereb

