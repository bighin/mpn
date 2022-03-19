#!/bin/bash
#
#SBATCH --job-name=mpn
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#

#
# ==== Info part ===== #
#
NOW=`date +%H:%M-%a-%d/%b/%Y`
echo '------------------------------------------------------'
echo 'This job is allocated on '$SLURM_JOB_CPUS_PER_NODE' cpu(s)'
echo 'Job is running on node(s): '
echo  $SLURM_JOB_NODELIST
echo '------------------------------------------------------'
echo 'WORKINFO:'
echo 'SLURM: job starting at           '$NOW
echo 'SLURM: sbatch is running on      '$SLURM_SUBMIT_HOST
echo 'SLURM: executing on cluster      '$SLURM_CLUSTER_NAME
echo 'SLURM: executing on partition    '$SLURM_JOB_PARTITION
echo 'SLURM: working directory is      '$SLURM_SUBMIT_DIR
echo 'SLURM: current home directory is '$(getent passwd $SLURM_JOB_ACCOUNT | cut -d: -f6)
echo ""
echo 'JOBINFO:'
echo 'SLURM: job identifier is         '$SLURM_JOBID
echo 'SLURM: job name is               '$SLURM_JOB_NAME
echo ""
echo 'NODEINFO:'
echo 'SLURM: number of nodes is        '$SLURM_JOB_NUM_NODES
echo 'SLURM: number of cpus/node is    '$SLURM_JOB_CPUS_PER_NODE
echo 'SLURM: number of cpus/task is    '$SLURM_CPUS_PER_TASK
echo 'SLURM: number of gpus/node is    '$SLURM_GPUS_PER_NODE
echo '------------------------------------------------------'
#
# ==== End of Info part ===== #
#

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

MPN_ADDITIONAL_ID=_0a ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0b ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0c ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0d ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0e ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0f ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0g ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0h ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0i ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0j ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0k ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0l ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0m ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0n ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0o ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0p ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0q ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0r ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0s ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0t ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0u ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0v ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0w ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0x ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0y ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_0z ./build/mpn $INIFILE &
wait
MPN_ADDITIONAL_ID=_1a ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1b ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1c ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1d ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1e ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1f ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1g ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1h ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1i ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1j ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1k ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1l ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1m ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1n ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1o ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1p ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1q ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1r ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1s ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1t ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1u ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1v ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1w ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1x ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1y ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_1z ./build/mpn $INIFILE &
wait
MPN_ADDITIONAL_ID=_2a ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2b ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2c ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2d ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2e ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2f ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2g ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2h ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2i ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2j ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2k ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2l ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2m ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2n ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2o ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2p ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2q ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2r ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2s ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2t ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2u ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2v ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2w ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2x ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2y ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_2z ./build/mpn $INIFILE &
wait
MPN_ADDITIONAL_ID=_3a ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3b ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3c ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3d ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3e ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3f ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3g ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3h ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3i ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3j ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3k ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3l ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3m ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3n ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3o ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3p ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3q ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3r ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3s ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3t ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3u ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3v ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3w ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3x ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3y ./build/mpn $INIFILE &
MPN_ADDITIONAL_ID=_3z ./build/mpn $INIFILE &
wait

NOW=`date +%H:%M-%a-%d/%b/%Y`
echo 'SLURM: job ending at             '$NOW
