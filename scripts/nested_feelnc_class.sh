#!/bin/bash --login
#SBATCH --account=pawsey1157-gpu
#SBATCH --nodes=1               # each subtask uses 1 node
#SBATCH --ntasks=1              # 1 subtask per file in the array-subtask
#SBATCH --time=24:00:00         # time per subtask
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --array=0-12

#############################################################
# Change to classify the clustered transcripts instead of the samples
# it should greatly reduce the resources needed
#############################################################

# Define directories
DIR="/scratch/pawsey1157/modanilevicz/sugarcane/nxf_work/data_stages"
SOFT_DIR="/software/projects/pawsey1157/modanilevicz/setonix/GitHub/sugarcane_lncRNA"
REF_ANNO="/scratch/pawsey1157/modanilevicz/sugarcane/reference/SofficinarumxspontaneumR570_771_v2.1.gene_exons.modified.gtf"

# Find project list
PROJECT_LIST=($(awk -F ',' '{gsub(/\r/, "", $5); print $5}' ${SOFT_DIR}/search/sugarcane_runs_cleaned.csv | sort -u | grep -v "Bioproject"))
PROJECT=${PROJECT_LIST[$SLURM_ARRAY_TASK_ID]}
echo "Processing $PROJECT"

# Load module and image
module load singularity/4.1.0-slurm
FEELnc="/software/projects/pawsey1157/groupResources/sharedImages/feelnc_0.1.1--r3.4.1_0.sif"

# 




	echo "Processing ${SAMPLE}"
	SAMPLE_DIR="${DIR}/candidate_lncRNA/${SAMPLE}"
	cd "${SAMPLE_DIR}" || exit

	split --lines=2048 "${SAMPLE_DIR}/lncRNAs.gtf" portioned_lncrna.
	
	NUMB_SPLITS=$(ls portioned_lncrna.* | wc -l)
	ARRAY_MAX=$((NUMB_SPLITS - 1))

set -exuo pipefail

PORTIONED_FILES=("${SAMPLE_DIR}"/portioned_lncrna.*)
FILE=\${PORTIONED_FILES[\$SLURM_ARRAY_TASK_ID]}

# Classify lncRNAs
srun -N 1 -n 1 --gres=gpu:4 \
    singularity exec \$FEELnc \
    FEELnc_classifier.pl \
    -i \${FILE} \
    -a \${REF_ANNO} \
    > \${FILE}_classified.txt

EOF

sbatch "${SOFT_DIR}/script/${SAMPLE}_array.sh"

COUNTER=$((COUNTER+1))
echo "Submitted sample ${SAMPLE}, ${COUNTER}/${TOTAL}"

done
