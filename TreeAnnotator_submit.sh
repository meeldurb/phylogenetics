#!/bin/bash

## man sbatch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=SLURM-LOGS/TreeAnnotator/TreeAnn-%A_%a.out

module load beast2/2.4.3


INPUT=$(ls $HOME/BEAST_results/*.Tree.trees.txt | awk 'NR=='$SLURM_ARRAY_TASK_ID'')
OUTDIR=$HOME/TreeAnnotator_results
OUTPUT=$(basename $INPUT).tre

cd $OUTDIR

treeannotator $INPUT $OUTDIR/$OUTPUT

