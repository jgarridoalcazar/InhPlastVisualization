#!/bin/bash
#Launching script for program.out
#Submission options:
#$-S /bin/bash
#$-N NESTSim
#$-cwd
#$-V

#./EDLUTKernel $BATCHARG
./EDLUTKernel -time $SIMTIME -nf $NETPREFIX$SGE_TASK_ID.dat -wf $WEIGHTFILE -if $INPUTFILE -log $OUTPUTPREFIX$SGE_TASK_ID.dat -ts $SIMULATIONTS
echo "Simulation Finished"
./SimilarityIndex $REFERENCEFILE $OUTPUTPREFIX$SGE_TASK_ID.dat $RESULTPREFIX$SGE_TASK_ID.dat $INITIAL_GRC $FINAL_GRC $INITIAL_TIME $INITIAL_TIME2 $BIN_SIZE $EVALUATION_TIME






#!/bin/bash
#$ -S /bin/bash
#
# set the P4_GLOBMEMSIZE
#$ -v P4_GLOBMEMSIZE=10000000
#
# Set the Parallel Environment and number of procs.
#$ -pe mpi 2

# Where we will make our temporary directory.
BASE="/tmp"

#
# make a temporary key
#
export KEYDIR=`mktemp -d $BASE/keys.XXXXXX`

#
# Make a temporary password.
# Makepasswd is quieter, and presumably more efficient.
# We must use the -s 0 flag to make sure the password contains no quotes.
#
if [ -x `which mkpasswd` ]; then
	export PASSWD=`mkpasswd -l 32 -s 0`
else
	export PASSWD=`dd if=/dev/urandom bs=512 count=100 | md5sum | gawk '{print $1}'`
fi

/usr/bin/ssh-keygen -t rsa1 -f $KEYDIR/tmpid -N "$PASSWD"

cat $KEYDIR/tmpid.pub >> $HOME/.ssh/authorized_keys

#
# make a script that will run under its own ssh-agent 
#
cat > $KEYDIR/launch-script <<"EOF"
#!/bin/bash
expect -c 'spawn /usr/bin/ssh-add $env(KEYDIR)/tmpid' -c \
	'expect "Enter passphrase for $env(LOGNAME)@$env(HOSTNAME)" \
		{ send "$env(PASSWD)\n" }' -c 'expect "Identity"'

echo

#
# Put your Job commands here.
#
#------------------------------------------------

/opt/openmpi/bin/mpirun -np $NSLOTS -machinefile $TMP/machines \
	/opt/hpl/openmpi-hpl/bin/xhpl

#/opt/mpich/gnu/bin/mpirun -np $NSLOTS -machinefile $TMP/machines \
	#/opt/hpl/mpich-hpl/bin/xhpl

#------------------------------------------------
EOF

chmod u+x $KEYDIR/launch-script

#
# start a new ssh-agent from scratch -- make it forget previous ssh-agent
# connections
#
unset SSH_AGENT_PID
unset SSH_AUTH_SOCK
/usr/bin/ssh-agent $KEYDIR/launch-script

#
# cleanup
#
grep -v "`cat $KEYDIR/tmpid.pub`" $HOME/.ssh/authorized_keys > $KEYDIR/authorized_keys
mv $KEYDIR/authorized_keys $HOME/.ssh/authorized_keys
chmod 644 $HOME/.ssh/authorized_keys

rm -rf $KEYDIR

