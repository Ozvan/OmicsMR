#!/bin/bash
#SBATCH --job-name=pwcoco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24GB
#SBATCH -t 12:00:00
#SBATCH -o ./log/output.%j.out
#SBATCH -e ./log/output.%j.err
#SBATCH --array=1-491%60
cd $SLURM_SUBMIT_DIR

list=$(head -n ${SLURM_ARRAY_TASK_ID} /scratch/richards/satoshi.yoshiji/11.pQTL/10.deCODE_cis_full_sumstats_batch/listbatch.txt | tail -n 1)
for protein in `cat $list`;
do
protein_name=`basename $protein ".cis.txt.gz"`  #in the format of "F11.2190_55"

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/bin/Rscript --vanilla 03.pwcoco.R $protein $protein_name
echo "done: $protein_name"
done
