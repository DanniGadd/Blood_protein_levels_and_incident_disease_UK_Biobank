#!/bin/bash -l
 
#SBATCH --job-name=batch
#SBATCH --cpus-per-task=5
#SBATCH --time 1-00:00:00
#SBATCH --mail-user=danni.gadd@biogen.com
#SBATCH --mail-type=END,FAIL

#SBATCH --array=1
module load R
srun Rscript 01_Protein_prep.R $SLURM_ARRAY_TASK_ID
