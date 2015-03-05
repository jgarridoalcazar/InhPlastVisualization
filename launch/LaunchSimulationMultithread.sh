#!/bin/sh
#$ -S /bin/sh
#$ -N sim-test-32
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -cwd
#$ -V
#$ -v OMP_NUM_THREADS=32
#$ -q 24Hbigmem
#$ -pe openmpi 32

python ./src/LaunchSimulation.py SimulationConfigTest.cfg

