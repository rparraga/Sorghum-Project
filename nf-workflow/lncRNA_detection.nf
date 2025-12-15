#! /usr/bin/env nextflow

/*
 * ========================================================================================
 * lncRNA Filtering Pipeline
 * ========================================================================================
 * This workflow filters a FASTA file of potential lncRNA candidates.
 *
 * Steps:
 * 1. Build a BLAST database from known rRNA and tRNA sequences.
 * 2. Remove sequences matching the rRNA/tRNA database using BLAST.
 * 3. Remove sequences matching known structured RNAs from the Rfam database using Infernal's cmscan.
 * 4. Keep only sequences predicted as non-coding by CPC2.
 *
 * Usage:
 * nextflow run <script_name>.nf \
 * --lncRNA_fasta /path/to/genome/reference \
 * --gtf /path/to/your/folder/of/candidates.gtf \
 * --outdir ./results
 * ========================================================================================
 */

// Define static parameters
params.genome = "something"
params.prjGTF = "/scratch/pawsey1157/modanilevicz/sugarcane/gtf"
params.prjDir = "/scratch/pawsey1157/modanilevicz/sugarcane/nxf_work/data_stages"


// Input channels
Channel
    .fromPath(params.sra_csv)
    .splitCsv(header : true)
    .map { row -> row["run_accession"]}
    .set {sample_id}


workflow {
    // Stage input files
    ch_lncRNA_fasta = file(params.lncRNA_fasta)
    ch_gtf = params.gtf ? file(params.gtf) : Channel.empty() // Create an empty channel if no GTF is provided

    // 1. Prepare reference database to remove tRNA and rRNA
    ch_blast_db = makeBlastDB(file(params.rRNA_tRNA_db))

    // 2. Blast against r/tRNA databases and filter them out
    ch_blast_filtered = blastFilter(ch_lncRNA_fasta, ch_gtf, ch_blast_db)

    // 3. Filter out known structured RNAs from Rfam
    ch_cmscan_filtered = cmscanFilter(ch_blast_filtered.fasta)

    // 4. Filter for non-coding potential using CPC2
    cpc2Filter(ch_cmscan_filtered)
}


// ##################################################
// ## PROCESSES
// ##################################################



process extractFASTA {
    input:
    val(sample_id)

    output:
    path("*fa")

    script:
    


}


process makeBlastDB {
    label 'makedb'
    tag "$sample_id"
    module "blast/2.12.0--pl5262h3289130_0"

    input:
    path rRNA_tRNA_sequences

    output:
    path "rRNA_tRNA.*"

    publishDir "${params.prjDir}/rRNA_tRNA_sequences", mode: 'copy'

    """
    makeblastdb -in ${rRNA_tRNA_sequences} -dbtype nucl -out rRNA_tRNA
    """
}

process blastFilter {
    label 'quality_checks'
    tag "Blast filtering ${sample_id}"
    module "blast/2.12.0--pl5262h3289130_0"

    input:
    path sample_fasta
    path db_files

    output:
    path "*.no_contaminants.fasta", emit: fasta
    path "*.no_rRNA.gtf", emit: gtf

    """
    blastn -query ${sample_fasta} -db ${db_files} -evalue 1e-5 -perc_identity 90 -outfmt 6 -out blast.out -num_threads ${task.cpus}
    
    # Get a list of sequence IDs that are contaminants
    cut -f1 blast.out | sort | uniq > contaminants.txt

    # Filter the FASTA file to remove contaminant sequences
    seqtk subseq -v ${lncRNA_fasta} contaminants.txt > lncRNAs_no_contaminants.fasta

    # Clear contaminants from the GTF
    grep -v -F -f contaminants.txt ${gtf} > ${lncRNA_fasta.baseName}_no_contaminants.gtf

    """
}

process cmscanFilter {
    input:
    path fasta from blastFilter.out

    output:
    path "lncRNAs.no_structured.fasta"

    """
    cmscan --cut_ga --rfam --nohmmonly --tblout hits.tblout ${params.rfam_cm} $fasta
    grep -v '^#' hits.tblout | awk '{print \$1}' | sort | uniq > structured_hits.txt
    seqtk subseq -v $fasta structured_hits.txt > lncRNAs.no_structured.fasta
    """
}

process cpc2Filter {
    input:
    path fasta from cmscanFilter.out

    output:
    path "lncRNAs.filtered.fasta"

    """
    cpc2.py -i $fasta -o cpc2_results.txt
    awk -F '\\t' '$8 == "noncoding" {print \$1}' cpc2_results.txt > noncoding_ids.txt
    seqtk subseq $fasta noncoding_ids.txt > lncRNAs.filtered.fasta
    """
}

workflow {
    makeBlastDB()
    blastFilter()
    cmscanFilter()
    cpc2Filter()
}
