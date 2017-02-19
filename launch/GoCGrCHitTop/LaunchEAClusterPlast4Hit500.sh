#!/bin/sh
#$ -S /bin/sh
#$ -N deap4-500-top-all-learning
#$ -m bea
#$ -M jesusgarrido@ugr.es
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -R y
#$ -cwd
#$ -v OMP_NUM_THREADS=1
#$ -q 72H
#$ -pe openmpi 128
#$ -hold_jid deap4-500-top-all-learning

export PATH=/SCRATCH/TIC117/jesusgarrido/autotools/bin:/SCRATCH/TIC117/jesusgarrido/NEST/nest210/insNoMPI/bin:/usr/local/apps/python-2.7.6/bin:$PATH
export PYTHONPATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest210/insNoMPI/lib/python2.7/site-packages/:/SCRATCH/TIC117/jesusgarrido/scipy/ins/lib/python2.7/site-packages/:$PYTHONPATH
export LD_LIBRARY_PATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest210/insNoMPI/lib/nest/:$LD_LIBRARY_PATH

mpirun-test -np $NSLOTS -ppn 8 python ./src/LaunchEvolutionaryAlgorithmMPI.py ./config/GoCGrCHitTop/EASearchGranularPlast4ClusterHit500.cfg
