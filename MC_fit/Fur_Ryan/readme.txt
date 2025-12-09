General cluster stuff
-module load scipy-stack (whenever you open a cluster window)
-sq (check progress of jobs)
-seff <job_ID> check how long it took, and why it ended (crashed completed etc). It tells you if the job didn't have enough memory or time etc, errors in python only show up in the slurm file...
-emacs slurm-<job_id> to see the output you get from a normal run including any bug reports that aren't related to messing up the cluster submission

"hatch.sh"
-how to run: ./hatch.sh <threshold> <one_sigma> example "./hatch.sh 80 10"
-you probably handle the threshold differently so may want to remove the two command line arguments from this script
-change the first line "jobName" variable to whatever you want the directory the results get saved in to be called
-calls the "job_fitnest.sh" file which submits the job...

"job_fitnest.sh"
-Don't run this file. It runs when hatch.sh runs it.
-"#SBATCH --time=7:59:00" probably the most important line, it's the time in hours:minutes:seconds you want to use. The longer the job the slower it takes to start, but you do not want the job to run out of time, so add 20% to how long you'd expect it to take, maybe +50% or more for the first time running on a cluster
-keep "#rm slurm-* || ls" commented out when debugging, but can comment it back in if storage space is a concern
-These two lines will probably need to be changed for you:
"cp fitPhoto5.py ${SUBDIR}/fitPhoto5.py
python3 fitPhoto5.py $1 $2  ${WORKDIR}/"
-fitPhoto5.py >> the python script that runs the calibration
-the code makes 50 numbered subdirectories which each have the results from a single toy dataset inside them. I do this using the ${WORKDIR}/ which gets saved as a prefix and whenever I save anything from a toy it is saved as ${WORKDIR}/<filename>

"fitPhoto5.py"
-This is my current version of the python script. Use your version, but maybe adapt it so it saved each output to a new spot with a command line ${WORKDIR}/ in "job_fitnest.sh"

"readWIMP.py/readCEvNS.py"
-Makes the style of plot with CEvNS or WIMP spectra in it.
-How to run: python readWIMP.py <output_directory> <threshold> <one_sigma>, example: "python readWIMP.py Quick_Test 160 20"
-This file is cursed. Do not try to figure out what numbers it regurgitates or you will lose your mind. I am sorry.
python readWIMP_v2.py /home/runze/Documents/results/TN_sims_2E5_pn_lead_low 400 50