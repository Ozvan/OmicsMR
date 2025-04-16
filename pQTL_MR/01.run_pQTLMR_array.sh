#!/bin/bash
#SBATCH --job-name=MR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24GB
#SBATCH -t 12:00:00
#SBATCH -o ./log/output.%j.out
#SBATCH -e ./log/output.%j.err
#SBATCH --array=1-189%50
cd $SLURM_SUBMIT_DIR

################
# This script submit array jobs
# Because many HPCs do not allow submitting thousands of jobs at one,
# we create a batch file each of which contain 10 proteins' path
# For example, batch_aa contains paths to A1BG.16561_9.tsv, /A4GALT.8759_29.tsv,... and ACP1.3858_5.tsv
# listbatch.txt is a list of these batch files
################

list=$(head -n ${SLURM_ARRAY_TASK_ID} /scratch/richards/satoshi.yoshiji/11.pQTL/03.cispQTL_decode2021_sep_batch/listbatch.txt| tail -n 1) 
for protein in `cat $list`;
do
protein_name=`basename $protein ".tsv"`  #in the format of "F11.2190_55"

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/bin/Rscript --vanilla 01.pQTL_MR.R $protein $protein_name
echo "done: $protein_name"
done
