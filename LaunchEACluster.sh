#!/bin/sh
#$ -S /bin/sh
#$ -N deap-7param2
#$ -m bea
#$ -M jesusgarrido@ugr.es
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -cwd
#$ -v OMP_NUM_THREADS=1
#$ -q 24H
#$ -pe openmpi 128
#$ -hold_jid deap-7param

export PATH=/SCRATCH/TIC117/jesusgarrido/autotools/bin:/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/bin:/usr/local/apps/python-2.7.6/bin:$PATH
export PYTHONPATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/lib/python2.7/site-packages/:/SCRATCH/TIC117/jesusgarrido/scipy/ins/lib/python2.7/site-packages/:$PYTHONPATH
export LD_LIBRARY_PATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/lib/nest/:$LD_LIBRARY_PATH

#python -m scoop --prolog /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/.env_var_ld -vvv -n 256 ./src/LaunchEvolutionaryAlgorithm.py SimulationConfigEACluster.cfg
mpirun -np 128 python ./src/LaunchEvolutionaryAlgorithmMPI.py SimulationConfigEACluster.cfg
