#!/bin/sh
#$ -S /bin/sh
#$ -N uego-7p8p8c20t
#$ -m bea
#$ -M jesusgarrido@ugr.es
#$ -o /SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer/results
#$ -j y
#$ -R y
#$ -cwd
#$ -v OMP_NUM_THREADS=1
#$ -q 72H
#$ -pe openmpi 128
#$ -hold_jid uego-7p8p8c20t

export PATH=/SCRATCH/TIC117/jesusgarrido/autotools/bin:/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/bin:/usr/local/apps/python-2.7.6/bin:$PATH
#export PYTHONPATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/lib/python2.7/site-packages/:/SCRATCH/TIC117/jesusgarrido/scipy/ins/lib/python2.7/site-packages/:/usr/local/apps/alhambra/python-2.7.6/lib/python2.7:$PYTHONPATH
export PYTHONPATH=/home/TIC117/jesusgarrido/.local/lib/python2.7/site-packages/readline-6.2.4.1-py2.7-linux-x86_64.egg:/home/TIC117/jesusgarrido/.local/lib/python2.7/site-packages/argparse-1.3.0-py2.7.egg:/home/TIC117/jesusgarrido/.local/lib/python2.7/site-packages/scoop-0.7.1.release-py2.7.egg:/usr/local/apps/alhambra/python-2.7.6/lib/python2.7/site-packages/setuptools-3.5.1-py2.7.egg:/usr/local/apps/alhambra/python-2.7.6/lib/python2.7/site-packages/pyqi-0.3.1-py2.7.egg:/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/lib/python2.7/site-packages:/SCRATCH/TIC117/jesusgarrido/scipy/ins/lib/python2.7/site-packages:/SCRATCH/TIC117/jesusgarrido/SpikingGranularLayer:/usr/local/apps/alhambra/python-2.7.6/lib/python27.zip:/usr/local/apps/alhambra/python-2.7.6/lib/python2.7:/usr/local/apps/alhambra/python-2.7.6/lib/python2.7/plat-linux2:/usr/local/apps/alhambra/python-2.7.6/lib/python2.7/lib-tk:/usr/local/apps/alhambra/python-2.7.6/lib/python2.7/lib-old:/usr/local/apps/alhambra/python-2.7.6/lib/python2.7/lib-dynload:/home/TIC117/jesusgarrido/.local/lib/python2.7/site-packages:/usr/local/apps/alhambra/python-2.7.6/lib/python2.7/site-packages
export LD_LIBRARY_PATH=/SCRATCH/TIC117/jesusgarrido/NEST/nest24/insNoMPI/lib/nest/:$LD_LIBRARY_PATH

mpirun -np $NSLOTS ./src/UEGO/Parallel/uego ./uegoconf/uegoini7p8p8c.optiParamCereb
