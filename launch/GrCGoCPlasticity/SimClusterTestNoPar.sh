#!/bin/sh
#$ -S /bin/sh
#$ -N TestNoParalela600
#$ -m bea
#$ -M jesusgarrido@ugr.es
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -R y
#$ -cwd
#$ -v OMP_NUM_THREADS=1
#$ -q NOParalela
# -pe openmpi 128
# -hold_jid deap4-500-top-no-learning

source $HOME/install/nest/ins/nompi/bin/nest_vars.sh

export PATH=$HOME/.local/bin:$HOME/install/cmake/ins/bin:$HOME/install/openssl/ins/bin:$HOME/install/cpython/ins/bin:$PATH
export LD_LIBRARY_PATH=$HOME/install/gsl/ins/lib:$HOME/install/openssl/ins/lib:$HOME/install/cpython/ins/lib:$HOME/install/readline/ins/lib:$HOME/install/libtool/ins/lib:$NEST_INSTALL_DIR/lib64:$LD_LIBRARY_PATH


python ./src/LaunchSimulation.py -c ./config/GrCGoCPlasticity/SimTestGrCGoC.cfg
