#!/usr/bin/env bash
# runs the monomer part of a DALTON PDE potential calculation
#
# options to be substituted
#  WORK_DIR: $WORK_DIR
#  SCRACTH : $SCRATCH
#  PROGPATH: $PROGPATH
#  JOB     : $JOB
#  MEMORY  : $MEMORY
#  NCPUS   : $NCPUS
#  LOPROP  : $LOPROP
#  MULMOM  : $MULMOM
#  POLMOM  : $POLMOM
#

export PATH=$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH

# compute the nescessary integrals through DALTON before we move on
export DALTON_TMPDIR=$SCRATCH
export DALTON_NUM_MPI_PROCS=$NCPUS
export OMP_NUM_THREADS=1

# in SLURM we prefer to use srun
#if [ -z ${SLURM_JOB_NAME+x} ]
#then
#    export DALTON_LAUNCHER="srun"
#fi

# run the calculation

# the .pot file contains coordinates of ALL polarizable sites
if [ -e $WORK_DIR/$JOB.h5 ]
then
    # if the .h5 file exists we will with a calculation using it
    if [ ! -e $WORK_DIR/$JOB.out ]
    then
        # but only if there is no output file
        $PROGPATH/dalton -mb $MEMORY -d -noarch -nobackup -noappend -get '$JOB.h5' -put '$JOB.h5' -o $JOB.out -dal $JOB.dal -pot temp.pot > $JOB.dalout

        if [ -e ${JOB}__temp.$JOB.h5 ]
        then
            mv ${JOB}__temp.$JOB.h5 $JOB.h5
        else
            echo "Output from PDE monomer calculation $JOB is missing."
        fi
    else
        echo "Skipping $JOB because output exists."
    fi
else
    echo "$JOB.h5 could not be found. Aborting."
fi

