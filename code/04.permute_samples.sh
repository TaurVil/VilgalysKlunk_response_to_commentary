#!/bin/bash
#SBATCH --get-user-env

index=${SLURM_ARRAY_TASK_ID}

module load R

sed -e s/ITERATION/${index}/g 04.section1_get_fst.R > tmp.$index.R
Rscript tmp.$index.R
rm tmp.$index.R


sed -e s/ITERATION/${index}/g 04.section2_get_pvals.R | sed -e s/TYPE/ML/g > tmp2.$index.b.R
Rscript tmp2.$index.R
rm tmp2.$index.R

sed -e s/ITERATION/${index}/g 04.section3_enrichment.R > tmp3.$index.R
Rscript tmp3.$index.R
rm tmp3.$index.R
echo DONE WITH $index
