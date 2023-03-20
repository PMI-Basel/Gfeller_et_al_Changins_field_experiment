#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --mem=40g
#SBATCH --output=run.out
#SBATCH --error=run.error
#SBATCH --job-name=ASV_clustering
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#load R module
module load foss/2018b
module load R/4.0.0-foss-2018b

#start scripts
taxa=$(awk '{print $2}' ../cmd/design.tab | sort | uniq)
bacteria=$(echo $taxa | grep -c b)
fungi=$(echo $taxa | grep -c f)
if [ ${bacteria} = 1 ]; then srun ./02.1_bacteria_ASV_clustering.R; fi
if [ ${fungi} = 1 ]; then srun ./02.2_fungi_ASV_clustering.R; fi
