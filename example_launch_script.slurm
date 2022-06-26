#!/bin/bash
#SBATCH --job-name=cycle1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nmovshov@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=64:00:00

. /home/nmovshov/miniconda3/etc/profile.d/conda.sh
conda activate

python drive-ppwd.py example_seeds/ppwd_seeds/neptune_N13_N2b.txt 4000 Neptune --toforder=4 -w=120 --ncores=30 --prefix=cycle1_ -j 2 4 >> rlog.$SLURM_JOB_ID.log