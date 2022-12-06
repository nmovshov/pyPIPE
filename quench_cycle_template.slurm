#!/bin/bash
#SBATCH --job-name=quench1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nmovshov@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --exclude=node001,node002,node003

. /home/nmovshov/miniconda3/etc/profile.d/conda.sh
conda activate

prefix=cycle1
rundir="$prefix"_saturn_run
seedfile=saturn/samples/deg4_tof4_J4_obs0_PF.txt
obs=Saturn1
dof=3
rtic=$SECONDS

for i in {1..1}
do
tic=$SECONDS
rm -rf $rundir

python drive-ppwd.py "$seedfile" 10 $obs -w=160 --ncores=40 -j 2 4 --toforder=4 --prefix="$prefix" -T=1e3 >> rlog.$SLURM_JOB_ID.log

python reseed.py "$rundir"/final.npz >> rlog.$SLURM_JOB_ID.log

python drive-ppwd.py fromreseed 10 $obs -w=160 --ncores=40 -j 2 4 --toforder=4 --prefix="$prefix" -T=1e2 >> rlog.$SLURM_JOB_ID.log

python reseed.py "$rundir"/final.npz >> rlog.$SLURM_JOB_ID.log

python drive-ppwd.py fromreseed 10 $obs -w=160 --ncores=40 -j 2 4 --toforder=4 --prefix="$prefix" -T=1e1 >> rlog.$SLURM_JOB_ID.log

python reseed.py "$rundir"/final.npz >> rlog.$SLURM_JOB_ID.log

python drive-ppwd.py fromreseed 100 $obs -w=160 --ncores=40 -j 2 4 --toforder=4 --prefix="$prefix" -T=1e0 >> rlog.$SLURM_JOB_ID.log

python save_ppwd_winners.py "$rundir"/final.npz wins_"$prefix".txt $dof $seedfile >> rlog.$SLURM_JOB_ID.log

toc=$SECONDS
echo "Elapsed time $((toc-tic)) seconds." >> rlog.$SLURM_JOB_ID.log

done
rtoc=$SECONDS
echo "Total runtime $(((rtoc-rtic)/3600)) hours." >> rlog.$SLURM_JOB_ID.log
