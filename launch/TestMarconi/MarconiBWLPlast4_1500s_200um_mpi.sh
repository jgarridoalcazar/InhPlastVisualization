#!/bin/bash
#PBS -l walltime=0:30:00
#PBS -l select=1:ncpus=1:mpiprocs=2:mem=3GB
#PBS -o /marconi_scratch/userexternal/jgarrido/SpikingGranularLayer/results/
#PBS -e /marconi_scratch/userexternal/jgarrido/SpikingGranularLayer/results/
#PBS -A Ppp27_3722
#PBS -N 200um-1th-2mpi
#PBS -q route

export PATH=/marconi/home/userexternal/jgarrido/autotools/bin:$PATH

cd $CINECA_SCRATCH
source ./nest/nest210/ins_mpi/bin/nest_vars.sh

module load intel
module load python/2.7.12
module load intelmpi/2017--binary
module load mkl/2017--binary
module load numpy/1.11.2--python--2.7.12
module load mpi4py/2.0.0--python--2.7.12
module load gsl/2.2.1--intel--pe-xe-2017--binary
module load scipy/0.18.1--python--2.7.12

cd SpikingGranularLayer
mpirun python ./src/LaunchSimulation.py -c ./config/TestMarconi/AllPlast4Hit1Th200um_mpi.cfg 2> ./results/$PBS_JOBNAME.$PBS_JOBID.txt > ./results/$PBS_JOBNAME.$PBS_JOBID.out
