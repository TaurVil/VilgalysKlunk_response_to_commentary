#!/bin/bash
#SBATCH --get-user-env

module load R
X=${SLURM_ARRAY_TASK_ID}

MYSITE=`head -$X sites_maf_to_permute.txt | tail -1`
sed -e s/MYSITE/$MYSITE/g 12a.get_1million_perms.R > g.$X.R

echo START $MYSITE

# run 100 iterations each of which will return 10k perms, for a total of 1 million
# avoids having unnecessarily large R scripts running, and would allow for additional permutations to be added easily
for f in `seq 1 100`; do echo START $f; Rscript ./g.$X.R; done 

bgzip perms.$MYSITE.txt

rm g.$X.R
echo $MYSITE
