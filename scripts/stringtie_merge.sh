#!/bin/bash --login

#SBATCH --account=pawsey1157
#SBATCH --output=arraystringtie-%j.out   # output of each task
#SBATCH --array=0-25           # match the number of input files
#SBATCH --nodes=1               # each subtask uses 1 node
#SBATCH --ntasks=1              # 1 subtask per file in the array-subtask
#SBATCH --cpus-per-task=32      # cpus per subtask
#SBATCH --mem=64GB              # Needed memory per array-subtask 
#SBATCH --time=6:00:00         # time per subtask

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
## Projects from file
SOURCE="/software/projects/pawsey1157/${USER}/setonix/GitHub/Sorghum-Project/search/sorghum_runs_projects_merged_filtered.csv"
PROJECT_LIST=($(awk -F ',' '{gsub(/\r/, "", $2); print $2}' $SOURCE | sort -u | grep -v "Bioproject"))
PROJECT_ID=${PROJECT_LIST[$SLURM_ARRAY_TASK_ID]}
echo "My input project is ${PROJECT_ID}"

## Generate a list of samples per project
PROJECT_FILE="/software/projects/pawsey1157/${USER}/setonix/GitHub/Sorghum-Project/search/${PROJECT_ID}.idx"
touch $PROJECT_FILE

SAMPLE_LIST=($(grep ${PROJECT_ID} $SOURCE | awk -F ',' '{print $1}'))

for sample in "${SAMPLE_LIST[@]}";
do echo "/scratch/pawsey1157/${USER}/sorghum/nxf_work/data_stages/stringtie_1stPass/${sample}.gtf" >> $PROJECT_FILE";
done

## Set environment variables
DIR="/scratch/pawsey1157/${USER}/sorghum/nxf_work/data_stages"
IMAGE="/software/projects/pawsey1157/groupResources/sharedImages/stringtie_3.0.0--h29c0135_0.sif"

module load singularity/4.1.0-slurm

#---
GENOME="/scratch/pawsey1157/rparraga/Sbicolor_730_v5.1.modified.gff3"
MERGED_GTF="${DIR}/stringtie_2Pass/${PROJECT_ID}.gtf"

# PROJECT STRINGTIE MERGE
srun -N 1 -n 1 -c 32 \
    singularity run $IMAGE \
    stringtie --merge -p 32 \
    -G "$GENOME" \
    -o "$MERGED_GTF" \
    "$PROJECT_FILE"

echo "FINISHED RUNNING"
