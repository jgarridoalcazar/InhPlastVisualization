#!/bin/sh
#$ -S /bin/sh
#$ -N deap-2param-2
#$ -hold_jid deap-2param
#$ -m bea
#$ -M jesusgarrido@ugr.es
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -cwd
#$ -V
#$ -v OMP_NUM_THREADS=1
#$ -q 24H
#$ -pe openmpi 256

python -m scoop --prolog /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/.env_var_ld -vvv -n 256 ./src/LaunchEvolutionaryAlgorithm.py SimulationConfigEACluster.cfg
