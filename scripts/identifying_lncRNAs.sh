#!/bin/bash --login
#SBATCH --account=pawsey1157
#SBATCH --output=arraylncID-%j.out   # output of each task
#SBATCH --array=0-418         	# match the number of input files
#SBATCH --nodes=1               # each subtask uses 1 node
#SBATCH --ntasks=1              # 1 subtask per file in the array-subtask
#SBATCH --time=20:00:00         # time per subtask
#SBATCH --cpus-per-task=32
#SBATCH --mem=100GB
#SBATCH --partition=work

set -exuo pipefail

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
FILE_LIST=($( awk -F ',' '{gsub(/\r/, "", $1); print $1}' /software/projects/pawsey1157/modanilevicz/setonix/GitHub/Sorghum-Project/search/sorghum_aligned_projects.csv | sort -u | grep -v "run_accession" ))
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID]}
echo ${FILE}

## Set environment variables
DIR="/scratch/pawsey1157/modanilevicz/sorghum/nxf_work/data_stages"
## Define the tool containers
CPC2="/software/projects/pawsey1157/groupResources/sharedImages/cpc2_latest.sif"
FEELnc="/software/projects/pawsey1157/groupResources/sharedImages/feelnc_0.1.1--r3.4.1_0.sif"
PLEK="/software/projects/pawsey1157/groupResources/sharedImages/plek_1.2--py310h8ea774a_10.sif"
GFFREAD="/software/projects/pawsey1157/groupResources/sharedImages/gffread_0.12.7--h077b44d_6.sif"

module load singularity/4.1.0-slurm
module load blast/2.12.0--pl5262h3289130_0

#---
## Set sample variables and run
GTF="${DIR}/stringtie_3Pass/${FILE}/${FILE}.gtf"

## Check reference annotation
REF_GENO="/scratch/pawsey1157/modanilevicz/sorghum/reference/Sbicolor_730_v5.0.fa"
REF_ANNO="/scratch/pawsey1157/modanilevicz/sorghum/reference/Sbicolor_730_v5.1.gene.gff3"
MOD_ANNO="${REF_ANNO/.gff3/.modified.gtf}"

if [[ ! -f "$MOD_ANNO" ]]; then
	srun -N 1 -n 1 -c 24 \
		singularity run $GFFREAD \
		gffread "$REF_ANNO" \
		-T --keep-genes -o "$MOD_ANNO" 

		sed -i 's/transcript_id /transcript_id "/g; s/; /"; /g; s/$/"/' "$MOD_ANNO"
		echo "Cleaned the ${MOD_ANNO}, all good to go."
		REF_ANNO="$MOD_ANNO"
else
	REF_ANNO="$MOD_ANNO"
	echo "Found the modified version of REF annotation, ${MOD_ANNO}, all good to go."
fi

mkdir -p "${DIR}/candidate_lncRNA/${FILE}/"
cd "${DIR}/candidate_lncRNA/${FILE}/" 

## FEELnc - Removing transcripts that overlap with mRNA exons and/or are shorten than 200 nt
echo "Start FEELnc"
srun -N 1 -n 1 -c 32 \
	singularity exec ${FEELnc} \
    	FEELnc_filter.pl \
    	-i $GTF \
    	-a $REF_ANNO \
	--monoex=-1 \
	--verbosity 0 \
	--proc=24 \
	> ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc.gtf
echo "Finished FEELnc"

# Verify the file was actually created and has content
if [ ! -f "${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc.gtf" ]; then
    echo "ERROR: FEELnc output is empty. Check if your input GTF and REF use the same Chromosome names."
    exit 1
fi

## Exclude rRNA & tRNAs
# Extract the fasta sequences
echo "Starting GFFread to generate fasta..."
srun -N 1 -n 1 -c 32 \
	singularity run $GFFREAD \
	gffread -w ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc.fa \
	-g ${REF_GENO} \
	${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc.gtf

# Make tRNA and rRNA databases
echo "Generating Blast database of t/rRNAs"
# Check that we have the files needed
if [[ ! -f "/scratch/pawsey1157/modanilevicz/sugarcane/contaminants/other_rnas.fa" ]]; then
	cd /scratch/pawsey1157/modanilevicz/sugarcane/contaminants/
	echo "Download the reference other RNA files"
	wget https://gtrnadb.org/genomes/eukaryota/Sbico301/sorBic3-tRNAs.fa
	wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz
	wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
	gunzip *gz
	cat *fa *fasta > other_rnas.fa
fi

if [[ ! -f "/scratch/pawsey1157/modanilevicz/sugarcane/contaminants/other_rna.n*" ]]; then
	cd /scratch/pawsey1157/modanilevicz/sugarcane/contaminants/
	echo "Generate the other_rna database"
	srun -N 1 -n 1 -c 32 \
	makeblastdb -in other_rnas.fa  -dbtype nucl -out other_rnas
fi

# Compare to the tRNA and rRNA database
echo "Removing contaminant RNA types with Blast"
cd "${DIR}/candidate_lncRNA/${FILE}/"
srun -N 1 -n 1 -c 32 \
	blastn -db /scratch/pawsey1157/modanilevicz/sugarcane/contaminants/other_rnas \
	-query ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc.fa \
	-out known_trRNA_hits.txt \
	-num_threads 32 \
	-max_target_seqs 1 \
	-max_hsps 1 \
	-perc_identity 95 \
	-evalue 1e-10 \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# Filter these out of the candidate lncRNA GTF
if [[ -s known_trRNA_hits.txt ]]; then
	cut -f1 known_trRNA_hits.txt | sort | uniq > contaminant_ids.txt
	grep -v -F -w -f contaminant_ids.txt \
		${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc.gtf \
       		> ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc_clean.gtf
else 
       cp ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc.gtf  ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc_clean.gtf	
fi

# Generate a fasta version of our clean GTF
srun -N 1 -n 1 -c 32 \
        singularity run $GFFREAD \
        gffread -w ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc_clean.fa \
        -g ${REF_GENO} \
        ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc_clean.gtf

## Calculate protein coding potential
# CPC2
srun -N 1 -n 1 -c 32 \
	singularity exec $CPC2 \
	python3 /CPC2_standalone-1.0.1/bin/CPC2.py \
	-i ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc_clean.fa \
	-o ${DIR}/candidate_lncRNA/${FILE}/cpc2.out	

# PLEK - throws a srun error, but plek runs correctly
srun -N 1 -n 1 -c 32 \
	singularity exec $PLEK \
	python /usr/local/bin/PLEK.py \
	-fasta ${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc_clean.fa \
	-out ${DIR}/candidate_lncRNA/${FILE}/plek.out \
	-thread 32 || [ $? -eq 1 ]

# Collect the output from CPC2 and PLEK to get consensus lncRNA candidates
## 1. Parse CPC2 (Column 8 is label, 'noncoding')
awk '$8 == "noncoding" {print $1}' ${DIR}/candidate_lncRNA/${FILE}/cpc2.out.txt \
	| sort > ${DIR}/candidate_lncRNA/${FILE}/ids_cpc2.txt

## 2. Parse PLEK (Lines start with >ID, look for 'Non-coding')
grep "Non-coding" ${DIR}/candidate_lncRNA/${FILE}/plek.out \
	| sed 's/>//g' \
       	| awk '{print $3}' \
	| sort > ${DIR}/candidate_lncRNA/${FILE}/ids_plek.txt

## 3. Find Intersection (Consensus)
## Comm finds lines common to sorted files. -12 suppresses unique lines, leaving only common ones.
comm -12 ${DIR}/candidate_lncRNA/${FILE}/ids_cpc2.txt \
        ${DIR}/candidate_lncRNA/${FILE}/ids_plek.txt \
	> ${DIR}/candidate_lncRNA/${FILE}/consensus_lncRNA_ids.txt

## 4. Filter GTF and FASTA based on Consensus IDs
# Filter GTF
grep -F -w \
	-f ${DIR}/candidate_lncRNA/${FILE}/consensus_lncRNA_ids.txt \
	${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc_clean.gtf \
	> ${DIR}/candidate_lncRNA/${FILE}/lncRNAs.gtf

## 5. Filter FASTA (Using seqkit, ensuring it is loaded)
seqkit grep -f ${DIR}/candidate_lncRNA/${FILE}/consensus_lncRNA_ids.txt \
	${DIR}/candidate_lncRNA/${FILE}/candidate_lncRNA_feelnc_clean.fa \
	> ${DIR}/candidate_lncRNA/${FILE}/lncRNAs.fa

### Classify lncRNAs
##srun -N 1 -n 1 -c 32 \
##    singularity exec $FEELnc \
##    FEELnc_classifier.pl \
##    -i ${DIR}/candidate_lncRNA/${FILE}/lncRNAs.gtf \
##    -a $REF_ANNO \
##    > ${DIR}/candidate_lncRNA/${FILE}/lncRNA_classes.txt

###################################
# Save the output on Acacia
module load rclone/1.68.1

cd ${DIR}/candidate_lncRNA/${FILE}/
zip -r ${FILE}_lncRNA.zip . -x "${FILE}_lncRNA.zip"

rclone copy ${FILE}_lncRNA.zip pawsey1157:sorghumlncrna/candidatelncrna/
