#!/bin/sh
#$ -S /bin/sh
#$ -N FT1000FanIn
#$ -m bea
#$ -M jesusgarrido@ugr.es
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -R y
#$ -cwd
#$ -v OMP_NUM_THREADS=1
#$ -q 72H
#$ -pe impi 128
#$ -hold_jid FT1000FanIn

source $HOME/install/nest/ins/nompi/bin/nest_vars.sh

module load alhambra/gcc-7.3.0

export PATH=$HOME/.local/bin:$HOME/install/cmake/ins/bin:$HOME/install/openssl/ins/bin:$HOME/install/cpython/ins/bin:$PATH
export LD_LIBRARY_PATH=$HOME/install/gsl/ins/lib:$HOME/install/openssl/ins/lib:$HOME/install/cpython/ins/lib:$HOME/install/readline/ins/lib:$HOME/install/libtool/ins/lib:$NEST_INSTALL_DIR/lib64:$LD_LIBRARY_PATH


mpirun -np $(($NSLOTS+1)) python ./src/LaunchEvolutionaryAlgorithmMPI.py ./config/GoCFanIn/EAMFGoCFanIn1000.cfg
