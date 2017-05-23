#!/usr/bin/env bash
# runs a DALTON LoProp calculation
#
# options to be substituted
#  WORK_DIR: $WORK_DIR
#  SCRACTH : $SCRATCH
#  PROGPATH: $PROGPATH
#  JOB     : $JOB
#  MEMORY  : $MEMORY
#  NCPUS   : $NCPUS
#

if [ ! -e $WORK_DIR/$JOB.loprop ]
then
    # if the .loprop log file does not exist we first check to see if
    # the

    # if the log file from DALTON exists we move on
    if [ ! -e $WORK_DIR/$JOB.out ]
    then

        # compute the nescessary integrals through DALTON before we move on
        export DALTON_TMPDIR=$SCRATCH
        export DALTON_NUM_MPI_PROCS=$NCPUS
        export OMP_NUM_THREADS=1
        #export DALTON_LAUNCHER="srun"

        # run the calculation
        $PROGPATH/dalton -mb $MEMORY -d -noarch -nobackup -noappend -get 'AOONEINT DALTON.BAS SIRIFC AOPROPER RSPVEC' -ow $JOB.dal > $JOB.dalout

        # make filenames appropriate for loprop scripts
        mv $JOB.AOONEINT AOONEINT
        mv $JOB.DALTON.BAS DALTON.BAS
        mv $JOB.SIRIFC SIRIFC
        mv $JOB.AOPROPER AOPROPER
        mv $JOB.RSPVEC RSPVEC

        tar cfj loprop_input.tar.bz2 RSPVEC AOPROPER SIRIFC DALTON.BAS AOONEINT
        rm -f RSPVEC AOPROPER SIRIFC DALTON.BAS AOONEINT

    fi

    # if we have the files we need then we compute the LoProp properties we need
    # but be sure to bomb with a meaningful error message to help the user
    # as much as possible.
    # TODO: Fixup the programs below to make it possible to supply path to
    #       the loprop executable
    if [ -e loprop_input.tar.bz2 ]
    then
        tar xfj loprop_input.tar.bz2

        if [ -e AOONEINT -a -e DALTON.BAS -a -e SIRIFC -a -e AOPROPER -a -e RSPVEC ]
        then
            export PYTHONPATH=/Users/css/Programs/loprop-1.0/lib/python2.7/site-packages:$PYTHONPATH
            ~/Programs/loprop-1.0/bin/loprop.py -a 1 -l 0 -t . --decimal 9 > $JOB.loprop

            rm -f RSPVEC AOPROPER SIRIFC DALTON.BAS AOONEINT

        else
            echo "Could not run LoProp because the following files are missing:"
            if [ ! -e AOONEINT ]; then echo "  - AOONEINT"; fi;
            if [ ! -e DALTON.BAS ]; then echo "  - DALTON.BAS"; fi;
            if [ ! -e SIRIFC ]; then echo "  - SIRIFC"; fi;
            if [ ! -e AOPROPER ]; then echo "  - AOPROPER"; fi;
            if [ ! -e RSPVEC ]; then echo "  - RSPVEC"; fi;
        fi
    else
        echo "'loprop_input.tar.bz2' could not be found. Did the DALTON calculation crash?."
    fi
fi
