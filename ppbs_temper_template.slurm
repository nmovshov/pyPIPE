#!/bin/bash
#SBATCH --job-name=temper1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nmovshov@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --exclude=node001,node002,node003,node004

. /home/nmovshov/miniconda3/etc/profile.d/conda.sh
conda activate

prefix=cycle1
rundir="$prefix"_uranus_run
seedfile=uranus/ppbs_one/PF_tof4_J0_obs1_monoseed.txt
obs=Uranus_b1
Js="2"
strictJs="2"
dof=2
fixrot=1

rtic=$SECONDS
for i in {1..1}
do
tic=$SECONDS
rm -rf $rundir

python drive-ppbs.py "$seedfile" 1000 $obs -w=160 --ncores=40 -j $Js --toforder=4 --prefix="$prefix" -T=3e3 --fix-rot=$fixrot >> rlog.$SLURM_JOB_ID.log

python reseed.py "$rundir"/final.npz >> rlog.$SLURM_JOB_ID.log

#python drive-ppbs.py fromreseed 400 $obs -w=160 --ncores=40 -j $Js --toforder=4 --prefix="$prefix" -T=3e2 --fix-rot=$fixrot >> rlog.$SLURM_JOB_ID.log

#python reseed.py "$rundir"/final.npz >> rlog.$SLURM_JOB_ID.log

#python drive-ppbs.py fromreseed 200 $obs -w=160 --ncores=40 -j $Js --toforder=4 --prefix="$prefix" -T=3e1 --fix-rot=$fixrot >> rlog.$SLURM_JOB_ID.log

#python reseed.py "$rundir"/final.npz >> rlog.$SLURM_JOB_ID.log

python drive-ppbs.py fromreseed 1000 $obs -w=160 --ncores=40 -j $Js --toforder=4 --prefix="$prefix" -T=1e0 --fix-rot=$fixrot >> rlog.$SLURM_JOB_ID.log

python save_ppbs_winners.py "$rundir"/final.npz wins_"$prefix".txt $dof --fix-rot=$fixrot --seedfile="$seedfile" --J-strict=$strictJs -o=$obs >> rlog.$SLURM_JOB_ID.log

toc=$SECONDS
echo "Elapsed time $((toc-tic)) seconds." >> rlog.$SLURM_JOB_ID.log

done
rtoc=$SECONDS
echo "Total runtime $(((rtoc-rtic)/3600)) hours." >> rlog.$SLURM_JOB_ID.log
