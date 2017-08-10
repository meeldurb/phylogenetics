#!/bin/bash

## man sbatch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=SLURM-LOGS/BEAST-%A_%a.out

module load beast2/2.4.3

INPUT=$(ls $HOME/xml_outfiles/*xml | awk 'NR=='$SLURM_ARRAY_TASK_ID'')  

cd $HOME/BEAST_results

beast $INPUT
