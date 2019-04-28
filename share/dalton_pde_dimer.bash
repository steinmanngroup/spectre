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

# 
if [ ! -e $WORK_DIR/$JOB.h5 ]
then
    # if the .h5 file does not exist we attempt to copy the monomer file to the dimer file
    NEW_FILE=${JOB}.h5
    OLD_FILE=${NEW_FILE%_*}_monomer.h5
    if [ ! -e $WORK_DIR/${OLD_FILE} ]
    then
        echo "The old .h5 file '${OLD_FILE}' was not found."
        exit 1
    else
        cp $WORK_DIR/${OLD_FILE} $WORK_DIR/$JOB.h5
    fi
fi

if [ -e $WORK_DIR/$JOB.h5 ]
then
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
