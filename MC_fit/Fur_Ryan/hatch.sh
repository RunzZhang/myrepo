nameJob="TN_sims_2E5"

SUBDIR="/home/runze/Documents/results/TN_sims_2E5"

if [ ! -d $SUBDIR ]; then
    mkdir $SUBDIR
fi

cp job_fitnest.sh ${SUBDIR}/job_fitnest.sh
cp MCMC.py ${SUBDIR}/fitPhoto5.py

cd $SUBDIR

#sbatch job_fitnest.sh $1 $2
sbatch job_fitnest.sh

sleep 1

echo submitted
