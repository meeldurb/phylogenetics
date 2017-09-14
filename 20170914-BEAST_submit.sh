#!/bin/bash

## man sbatch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=SLURM-LOGS/BEAST-%A_%a.out

module load beast2/2.4.3

INPUT=$(ls $HOME/results_BEAUTi_xml/20170914-2duppairs_ntclans/*xml | awk 'NR=='$SLURM_ARRAY_TASK_ID'')  

cd $HOME/results_BEAST/20170914-BEAST_clans_ntseq

beast $INPUT
