#!/bin/bash
#
#PBS -l ncpus=320
#PBS -l walltime=00:30:00
#PBS -l mem=1280GB
#
#PBS -q normal
#PBS -P jh2
#
#PBS -M egentry@ucsc.edu
#PBS -m e
#
#PBS -e ./logs/
#PBS -o ./logs/
#
#PBS -S /bin/bash
#
#PBS -N 300_HD
#

cd $PBS_O_WORKDIR

source $HOME/.bashrc

# number of processors
NP=$(wc -l $PBS_NODEFILE | awk '{print $1}')

RUN_NAME=cluster_cooling_300
RUN_DIR=./runs/$RUN_NAME

INPUTS_DIR=$RUN_DIR/inputs

OUTPUT_FILE=./logs/${PBS_JOBNAME}.byhand.o${PBS_JOBID}

touch $OUTPUT_FILE

echo $PBS_NODEFILE >> $OUTPUT_FILE
cat $PBS_NODEFILE >> $OUTPUT_FILE

COUNTER=0
while [ ! -f $INPUTS_DIR/end ]; do
	echo Counter: $COUNTER 1>>$OUTPUT_FILE 2>&1
	let COUNTER+=1

	# default behavior: loop will break before the simulation gets started
	#   The simulation is only allowed to run if `prepare_for_restart.py` decides it's okay, and deletes the `end` file
	touch $INPUTS_DIR/end

	cd scripts
	echo "about to enter 'prepare_for_restart.py'" 1>>../$OUTPUT_FILE 2>&1
	python3 prepare_for_restart.py ../$INPUTS_DIR/end ../$RUN_DIR 1>>../$OUTPUT_FILE 2>&1
	echo "done with 'prepare_for_restart.py'" 1>>../$OUTPUT_FILE 2>&1
	cd $PBS_O_WORKDIR

	# make sure we don't run the code if the python failed, or decided the code shouldn't be run
	if [ -f $INPUTS_DIR/end ]; then
		echo "Found 'inputs/end' file -- breaking early" 1>>$OUTPUT_FILE 2>&1
		break
	fi

	RESTART_TYPE=2 # restart from snapshot file
	
	time mpirun -hostfile $PBS_NODEFILE -np $NP $INPUTS_DIR/GIZMO $INPUTS_DIR/$RUN_NAME.params.restart $RESTART_TYPE  1>>$OUTPUT_FILE 2>&1
	
done

###compress all resulting files, if the simulation ended successfully
# cd $RUN_DIR/outputs
# qsub $PBS_O_WORKDIR/scripts/compress_and_delete_all

# cd $PBS_O_WORKDIR

