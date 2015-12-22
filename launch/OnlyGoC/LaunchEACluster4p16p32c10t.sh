#!/bin/sh
#$ -S /bin/sh
#$ -N deap-4p16p32c10t
#$ -m bea
#$ -M jesusgarrido@ugr.es
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -R y
#$ -cwd
#$ -v OMP_NUM_THREADS=1
#$ -q 72H
#$ -pe openmpi 128
#$ -hold_jid deap-4p16p32c10t

export PATH=/SCRATCH/TIC117/jesusgarrido/autotools/bin:/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/bin:/usr/local/apps/python-2.7.6/bin:$PATH
export PYTHONPATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/lib/python2.7/site-packages/:/SCRATCH/TIC117/jesusgarrido/scipy/ins/lib/python2.7/site-packages/:$PYTHONPATH
export LD_LIBRARY_PATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/lib/nest/:$LD_LIBRARY_PATH

#python -m scoop --prolog /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/.env_var_ld -vvv -n 256 ./src/LaunchEvolutionaryAlgorithm.py SimulationConfigEACluster.cfg
mpirun -np $NSLOTS python ./src/LaunchEvolutionaryAlgorithmMPI.py SimulationConfigEACluster4p16p32c10t.cfg
