#!/bin/bash
#SBATCH --job-name=Sseed
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nmovshov@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=24:00:00

. /home/nmovshov/miniconda3/etc/profile.d/conda.sh
conda activate

python drive-ppwd.py saturn/seeds/deg4_quadratic_base_seed.txt 22000 Saturn -w=160 --ncores=40 -j=0 -s --prefix=period_flex --preserve-period=0 >> rlog.$SLURM_JOB_ID.log

