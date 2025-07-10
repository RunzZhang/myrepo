#!/bin/bash                                                                                                                                                    
#SBATCH --time=22:59:00
#SBATCH --array=1-50
#SBATCH --account=def-kenclark-ab                                                                                                                                             
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB

#module load python

#module purge
#module load StdEnv/2023
module load scipy-stack

#rm slurm-* || ls                                                                                                                                                                                         

WORKDIR=${SLURM_ARRAY_TASK_ID}

if [ ! -d $WORKDIR ]; then
    mkdir $WORKDIR
fi

#python3 fitPhoto5.py $1 $2  ${WORKDIR}/
python3 MCMC.py ${WORKDIR}/
#rm slurm-* || ls                                                                                                                                                                                           

#echo information

sleep 1
