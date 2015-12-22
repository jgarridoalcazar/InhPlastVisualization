#!/bin/sh
#$ -S /bin/sh
#$ -N sim-test-256
#$ -m bea
#$ -M jesusgarrido@ugr.es
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -cwd
#$ -V
#$ -v OMP_NUM_THREADS=16
#$ -q 72H
#$ -pe openmpi 256

export PATH=/SCRATCH/TIC117/jesusgarrido/autotools/bin:/SCRATCH/TIC117/jesusgarrido/NEST/nest28/insMPI/bin:/usr/local/apps/python-2.7.6/bin:$PATH
export PYTHONPATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest28/insMPI/lib/python2.7/site-packages/:/SCRATCH/TIC117/jesusgarrido/scipy/ins/lib/python2.7/site-packages/:$PYTHONPATH
export LD_LIBRARY_PATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest28/insMPI/lib/nest/:$LD_LIBRARY_PATH

mpirun-test -np 16 -ppn 1 python /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/src/LaunchSimulation.py -c /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/config/GoCGrC/SimulationConfigGranularTestCluster.cfg
#python /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/src/LaunchSimulation.py -c /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/config/GoCGrC/SimulationConfigGranularTestCluster.cfg

