#!/bin/bash
# pe request
#$ -pe mpi_4 4
#$ -N  pow
#$ -S /bin/bash
#$ -q all.q
#$ -V
#$ -cwd

echo "Got $NSLOTS slots."

    #######################################################
    ### openmpi 1.6.4 (w/ Intel compiler)
    ### >>> THIS IS THE LABORATORY STANDARD!
    ### >>>  IT IS NOT RECOMMENDED TO MODIFY
    ### >>>  THE DEFAULT SETTINGS AS SHOWN HERE.
    #
    export PSM_SHAREDCONTEXTS_MAX=`echo $PE | awk -F_ '{print $2}'`
    export PSM_RANKS_PER_CONTEXT=4 
    #
    #######################################################

    mpirun -n $NSLOTS /home/amir/.exe/seed_pow
