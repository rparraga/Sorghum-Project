#! /usr/bin/env nextflow
params.sra_csv = '/software/projects/pawsey1157/rparraga/setonix/GitHub/Sorghum-Project/Sorghum_Cleaned.csv'
params.prjDir = '/scratch/pawsey1157/rparraga/sorghum/nxf_work/data_stages' 
params.threads = 16

Channel
    .fromPath(params.sra_csv)
    .splitCsv(header : true)
    .collect() 
    .flatMap { rows ->
	def total_lines = rows.size()
	def valid_rows = rows.findAll { row ->
            def run_id = row.Run?.trim()
            run_id && (run_id.startsWith("SRR") || run_id.startsWith("ERR"))
	}	
	if (valid_rows.size() != total_lines) {
            error "Input file is invalid. Expected ${total_lines} SRR/ERR IDs, but found only ${valid_rows.size()} valid entries."
        }
	return valid_rows
     } 
    .map { row ->
        def sra_id = row.Run
        def is_paired = row.LibraryLayout.toLowerCase() == 'paired'
        tuple(sra_id, is_paired)
    }
    .set { sra_ids }

Channel
    .fromPath(params.sra_csv)
    .splitCsv(header : true) 
    .map {row ->
	def bioproject = row.BioProject
	}
    .unique()
    .set { bioproject_list }


/*
 * ================================================================================================
 * WORKFLOW
 * ================================================================================================
 */


workflow {
    downloadandConvertSRA(sra_ids)
//    fastpTrim(downloadandConvertSRA.out.reads)
//    fastqc(fastpTrim.out.trimmed_reads)
//    checkFastqcQuality(fastqc.out.fastqc_zips)
//
//    checkFastqcQuality.out.qc_summary
//        .collectFile(name: 'qc_summary.txt', storeDir: "${params.prjDir}/quality_control")
//
//    genome_index_ch = Channel.fromPath("${params.genome_library}*ht2l")
//                             .collect() 
//
//    align_out = align2genome(genome_index_ch, fastpTrim.out.trimmed_reads)
//    samtools(align_out.sam)
//    stringtie(samtools.out.bam)
//	
//    all_gtf = stringtie.out.gtf
//		.collect()			
//	 
//    stringtieMerge(bioproject_list, all_gtf)
//    stringtieQuantify()
//    getFasta()
}

/*
 * ================================================================================================
 * PROCESSES
 * ================================================================================================
 */

process downloadandConvertSRA {
    tag "$sra_id"
    publishDir "${params.prjDir}/fastq", mode: 'copy'

    input:
    tuple val(sra_id), val(is_paired)

    output:
    tuple val(sra_id), path("${sra_id}_*.fastq"), path("actual_is_paired.txt"), emit: reads

    script:
    def dump_cmd = is_paired ?
        "fasterq-dump $sra_id --split-files --outdir . --threads $params.threads" :
        "fasterq-dump $sra_id --outdir . --threads $params.threads"
     
    """
    set -euo pipefail

    echo "Downloading $sra_id (paired=${is_paired})"

    $dump_cmd
    
    NUM_FASTQ=\$(ls -1 *.fastq | wc -l)

    if [ "\$NUM_FASTQ" -eq 2 ]; then
	echo "true" > actual_is_paired.txt
    elif [ "\$NUM_FASTQ" -eq 1 ]; then
        echo "false" > actual_is_paired.txt
        mv *.fastq ${sra_id}_1.fastq
    else
        echo "Error: Expected 1 or 2 FASTQ files for $sra_id, but found \$NUM_FASTQ."
        exit 1
    fi
    """
}


//process fastpTrim {
//    label 'qualitychecks'
//    tag "$sample_id"
//    publishDir "${params.prjDir}/trimmed", mode: 'copy'
//
//    input:
//    tuple val(sample_id), path(reads), path(is_paired)
//    
//    output:
//    tuple val(sample_id), path("*_trimmed.fastq"), path(is_paired), emit:trimmed_reads
//    path("${sample_id}_fastp.html"), emit: html
//    path("${sample_id}_fastp.json"), emit: json
//
//    script:
//    def R1 = reads.find { it.name.contains('_1.fastq') }
//    def R2 = reads.find { it.name.contains('_2.fastq') }
//    def out_R1 = "${sample_id}_1_trimmed.fastq"
//    def out_R2 = "${sample_id}_2_trimmed.fastq"
//    
//    """
//    IS_PAIRED=\$(cat $is_paired)
//    echo "Processing $sample_id (paired=\${IS_PAIRED})"
//
//    if [ "\$IS_PAIRED" = "true" ]; then
//        fastp -i $R1 -I $R2 \
//              -o $out_R1 \
//              -O $out_R2 \
//              --detect_adapter_for_pe \
//              --html ${sample_id}_fastp.html \
//              --json ${sample_id}_fastp.json \
//              --thread $params.threads
//      
//    else 
//        fastp -i $R1 \
//              -o $out_R1 \
//              --html ${sample_id}_fastp.html \
//              --json ${sample_id}_fastp.json \
//              --thread $params.threads
//    fi
//    """
//}
//
//
//process fastqc {
//    tag "$sample_id"
//    module "fastqc/0.11.9--hdfd78af_1"
//    publishDir "${params.prjDir}/fastqc", mode: 'copy'
//
//    input:
//    tuple val(sample_id), path(trimmed_reads), path(is_paired)
//    
//    output:
//    tuple val(sample_id), path("*_fastqc.zip"), emit: fastqc_zips
//
//    script:
//    """
//    fastqc $trimmed_reads --threads $params.threads
//    """
//}
//
//
//process checkFastqcQuality {
//    tag "$sample_id"
//    label 'quality_checks'
//    publishDir "${params.prjDir}/fastqc_failed", mode: 'copy', pattern: "*failed_summary.txt"
//
//    input:
//    tuple val(sample_id), path(fastqc_zips)
//
//    output:
//    path "${sample_id}_qc_result.txt", emit: qc_summary
//
//
//    script:
//    """
//    mkdir ${sample_id}_qc_check
//
//    for zip in *_fastqc.zip; do
//        unzip "\$zip" -d ${sample_id}_qc_check
//    done
//
//    summary_files=\$(find ${sample_id}_qc_check -name summary.txt)
//
//    grep -E 'FAIL' \$summary_files > ${sample_id}_qc_check/failed_modules.txt || true
//    n_issues=\$(wc -l < ${sample_id}_qc_check/failed_modules.txt)
//
//    if [ "\$n_issues" -gt 3 ]; then
//        echo -e "$sample_id FAILED\\n----------------" > ${sample_id}_failed_summary.txt
//        cat ${sample_id}_qc_check/failed_modules.txt >> ${sample_id}_failed_summary.txt
//        echo -e "$sample_id\tFAIL"
//    else
//        echo -e "$sample_id\tPASS"
//    fi > ${sample_id}_qc_result.txt
//    """
//}
//
//
//process align2genome {
//    tag "$sample_id"
//    label 'genome_align'
//    publishDir "${params.prjDir}/hisat2_align", mode: 'copy'
//
//    input:
//    path(index_files)
//    tuple val(sample_id), path(trimmed_reads), path(is_paired)
//
//    output:
//    tuple val(sample_id), path("${sample_id}.sam"), emit: sam
//    path("${sample_id}_align.log"), emit: log
//    
//    script:
//    def genome_basename = "SofficinarumxspontaneumR570_771_v2.0"
//    def R1 = trimmed_reads.find { it.name.contains('_1_trimmed.fastq') }
//    def R2 = trimmed_reads.find { it.name.contains('_2_trimmed.fastq') }
//    
//    """	
//    IS_PAIRED=\$(cat $is_paired)
//    echo "Processing $sample_id (paired=\$IS_PAIRED)"
//    
//    if [ "\$IS_PAIRED" = "true" ]; then 
//    echo "running paired command"
//    hisat2 -p ${params.threads} \
//           --phred33 --no-unal \
//           -x ${genome_basename} \
//           -1 ${R1} -2 ${R2} \
//           -S ${sample_id}.sam \
//           2> ${sample_id}_align.log
//    
//    else
//    echo "running single command"
//    hisat2 -p ${params.threads} \
//           -x ${genome_basename} \
//           -U ${R1} \
//           -S ${sample_id}.sam \
//           2> ${sample_id}_align.log
//    fi  
//    """
//}
//
//process samtools {
//    tag "$sample_id"
//    publishDir "${params.prjDir}/samtools", mode: 'copy'
//
//    input:
//    tuple val(sample_id), path(sam)
//
//    output:
//    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
//
//    script:	
//    """
//    samtools sort -@ $params.threads \
//		-o ${sample_id}.bam \
//  		$sam
//    """
//}
//
//process stringtie {
//    tag "$sample_id"
//    label 'read_quant'
//    publishDir "${params.prjDir}/stringtie_1stPass", mode: 'copy'
//
//    input:
//    tuple val(sample_id), path(bam)
//
//    output:
//    tuple val(sample_id), path("${sample_id}.gtf"), emit: gtf
//
//    script:
//    """
//    stringtie  -p ${task.cpus} \
//                -G ${params.genome_gtf} \
//                -o ${sample_id}.gtf \
//                -l $sample_id \
//                $bam
//    """
//}
//
 
//process stringtieMerge {
//    tag "$sample_id"
//    label 'read_quant'
//    publishDir "${params.prjDir}/stringtie_merge", mode: 'copy'
//
//    input:
//    // tuple val(bioproject), path(bioproject_list) //file with sample gtf files generated earlier
//    // genome reference GFF
//
//    output:
//    tuple val(bioproject), path("${bioproject}.gtf", emit: merged
//
//    script:
//
//    """
//    // generate the bioproject file on the go
//
//    stringtie --merge -p ${task.cpus} \
//                -G "${params.genome_gtf}" \
//                -o "${bioproject}_merged.gtf" \
//                -l "\$FILE_LIST"
//    """
//}






























