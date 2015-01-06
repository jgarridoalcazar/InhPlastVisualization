#!/bin/sh
#$ -S /bin/sh
#$ -N deap-2param
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -cwd
#$ -V
#$ -v OMP_NUM_THREADS=1
#$ -q 24H
#$ -pe openmpi 128

python -m scoop --prolog /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/.env_var_ld -vvv -n 128 ./src/LaunchEvolutionaryAlgorithm.py SimulationConfigEACluster.cfg