nameJob="Fit5_Quick2_320"

SUBDIR=${PWD}/${nameJob}

if [ ! -d $SUBDIR ]; then
    mkdir $SUBDIR
fi

cp job_fitnest.sh ${SUBDIR}/job_fitnest.sh
cp fitPhoto5.py ${SUBDIR}/fitPhoto5.py

cd $SUBDIR

sbatch job_fitnest.sh $1 $2

sleep 1

echo submitted
