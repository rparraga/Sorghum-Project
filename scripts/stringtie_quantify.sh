#!/bin/bash --login
#SBATCH --account=pawsey1157
#SBATCH --output=arraystringtie-%j.out   # output of each task
#SBATCH --array=0-259          # match the number of input files
#SBATCH --nodes=1               # each subtask uses 1 node
#SBATCH --ntasks=1              # 1 subtask per file in the array-subtask
#SBATCH --cpus-per-task=32      # cpus per subtask
#SBATCH --mem=64GB              # Needed memory per array-subtask 
#SBATCH --time=4:00:00         # time per subtask

set -euo pipefail

#--- 
echo "All jobs in this array have:"
echo "- SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
echo "- SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX}"
  
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
 
#--- 
## Sample list
FILE_LIST=($(awk -F ',' '{gsub(/\r/, "", $4); print $4}' /software/projects/pawsey1157/modanilevicz/setonix/GitHub/sugarcane_lncRNA/search/sugarcane_runs_cleaned.csv | sort -u | grep -v "run_accession"))
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID]}
echo ${FILE}

## Find out the sample project
PROJECT=$(grep $FILE /software/projects/pawsey1157/modanilevicz/setonix/GitHub/sugarcane_lncRNA/search/sugarcane_runs_cleaned.csv | awk -F ',' '{gsub(/\r/, "", $5); print $5}')

## Set environment variables
DIR="/scratch/pawsey1157/modanilevicz/sugarcane/nxf_work/data_stages"
IMAGE="/software/projects/pawsey1157/groupResources/sharedImages/stringtie_3.0.0--h29c0135_0.sif"

module load singularity/4.1.0-slurm

#---
## Set sample variables and run
REF_GTF="${DIR}/stringtie_2Pass/${PROJECT}.gtf"

mkdir -p "${DIR}/stringtie_3Pass/${FILE}/"
OUT_GTF="${DIR}/stringtie_3Pass/${FILE}/${FILE}.gtf"

BAM="${DIR}/samtools/${FILE}.bam"

srun -N 1 -n 1 -c 32 \
    singularity run $IMAGE \
    stringtie -e -B -p 32 \
    -G  "${REF_GTF}" \
    -o  "${OUT_GTF}" \
    $BAM

