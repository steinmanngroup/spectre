#!/usr/bin/env bash
# runs a DALTON calculation
#
# options to be substituted
#  WORK_DIR: $WORK_DIR
#  SCRACTH : $SCRATCH
#  PATH    : $PROGPATH
#  JOB     : $JOB
#  MEMORY  : $MEMORY
#  NCPUS   : $NCPUS
#  PDEOPT  : $PDEOPT
#

export PATH=$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH

if [ ! -e $WORK_DIR/$JOB.out ]
then

    export DALTON_TMPDIR=$SCRATCH
    export DALTON_NUM_MPI_PROCS=$NCPUS
    export OMP_NUM_THREADS=1
    # in SLURM we prefer to use srun
    #if [ -z ${SLURM_JOB_NAME+x} ]
    #then
    #    export DALTON_LAUNCHER="srun"
    #fi
    
    # run the calculation
    $PROGPATH/dalton -mb $MEMORY $PDEOPT -d -noappend -noarch -o $JOB.out -dal $JOB.dal -pot $JOB.pot > $JOB.dalout

else
    echo "Skipping $JOB because output exists."
fi

