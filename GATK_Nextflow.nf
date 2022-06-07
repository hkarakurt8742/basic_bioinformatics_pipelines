#! /usr/bin/env nextflow

/* 
 * pipeline input parameters 
 */
params.reads = ""
params.genome_file = ""
params.outdir = ""
params.known_sites = ""
params.target_bed = ""

log.info """\
         ALIGNMENT AND VARIANT CALLING PIPELINE FOR LOCAL USE (Version 1)  
         Written by Hamza Umut Karakurt  
         ===================================
         Genome: ${params.genome_file}
         Reads        : ${params.reads}
         Outdir       : ${params.outdir}
         Known Sites        : ${params.known_sites}
         Target Bed       : ${params.target_bed}
         """
         .stripIndent()

 
/* 
 * define the `index` process that create a binary index 
 * given the transcriptome file
 */


read_pairs_ch = Channel .fromFilePairs(params.reads)
read_pairs2_ch = Channel .fromFilePairs(params.reads)


process alignment {
    cpus 14
    tag "BWA on $sample_id"
    publishDir "${params.outdir}", mode: 'copy'
         
    input:
    path genome from params.genome_file
    tuple val(sample_id), path(reads) from read_pairs_ch
 
    output:
    path "${sample_id}.bam" into aligned_ch, aligned2_ch
 
    script:
    """
    bwa mem -M -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:Illumina' -t $task.cpus ${params.genome_file} ${reads[0]} ${reads[1]} | samtools sort -@ $task.cpus -o ${sample_id}.bam -
    """
}


process bamprocess {
    cpus 14
    tag "$Indexing on sample_id"
    publishDir "${params.outdir}", mode: 'copy'
         
    input:
    file sample_id from aligned_ch
 
    output:
    path "${sample_id.baseName}.bam.bai" into bam_indexed_ch
    path "${sample_id.baseName}_stats.txt" into bam_stats_ch
 
    script:
    """
    samtools index -@ $task.cpus ${sample_id} ${sample_id.baseName}.bam.bai
    samtools stats ${sample_id} > ${sample_id.baseName}_stats.txt
    """
}

process fastqc {
    tag "FASTQC on $sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) from read_pairs2_ch

    output:
    path "fastqc_${sample_id}_logs" into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}  


process markduplicates {
    tag "MarkDuplicates on $sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sample_id from aligned2_ch

    output:
    file "${sample_id.baseName}.marked.bam" into marked_bams_ch, marked_bams2_ch, marked_bams3_ch
    file "${sample_id.baseName}_metrics.txt" into picard_metrics_ch
    


    script:
    """
    picard MarkDuplicates INPUT=${sample_id} OUTPUT=${sample_id.baseName}.marked.bam METRICS_FILE=${sample_id.baseName}_metrics.txt CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT
    """  
}  


process baserecalibrator {
    tag "BaseRecalibrator on $sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sample_id from marked_bams_ch
    path genome from params.genome_file
    path known_sites from params.known_sites

    output:
    file "${sample_id.baseName}_recal_data.table" into recal_tables_ch
    

    script:
    """
    /home/huk/Desktop/Tool/gatk-4.2.0.0/gatk BaseRecalibrator -I ${sample_id}  -R ${params.genome_file} --known-sites ${params.known_sites} -O ${sample_id.baseName}_recal_data.table
    """
}



process applybqsr {
    tag "ApplyBQSR on $sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sample_id from marked_bams2_ch
    path genome from params.genome_file
    path recal_table from recal_tables_ch

    output:
    file "${sample_id.baseName}.marked.recal.bam" into marked_recalled_bams_ch
    

    script:
    """
    /home/huk/Desktop/Tool/gatk-4.2.0.0/gatk ApplyBQSR -R ${params.genome_file} -I ${sample_id} --bqsr-recal-file ${recal_table} -O ${sample_id.baseName}.marked.recal.bam
    """
}


process haplotypecaller {
    cpus 14
    tag "HaploTypeCaller on $sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sample_id from marked_bams3_ch
    path genome from params.genome_file
    path known_sites from params.known_sites
    path target_bed from params.target_bed

    output:
    file "${sample_id.baseName}.vcf.gz" into vcf_outputs_ch
    

    script:
    """
    /home/huk/Desktop/Tool/gatk-4.2.0.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.genome_file} -I ${sample_id} -O ${sample_id.baseName}.vcf.gz --dbsnp ${params.known_sites} -L ${params.target_bed} --native-pair-hmm-threads $task.cpus
    """
}



process variantfilter {
    cpus 14
    tag "VariantFilter on $sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sample_id from vcf_outputs_ch
    path genome from params.genome_file

    output:
    file "${sample_id.baseName}_filtered.vcf.gz" into filtered_vcf_ch
    file "${sample_id}.tbi" into vcf_index
    

    script:
    """
       /home/huk/Desktop/Tool/gatk-4.2.0.0/gatk IndexFeatureFile -I ${sample_id}
    /home/huk/Desktop/Tool/gatk-4.2.0.0/gatk VariantFiltration -R ${params.genome_file} -V ${sample_id} -O ${sample_id.baseName}_filtered.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNPFilter"
    """
}



process annotation {
    conda "/home/huk/anaconda3/envs/vep"
    cpus 14
    tag "Annotation on $sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sample_id from filtered_vcf_ch

    output:
    file "${sample_id.baseName}_vep.txt" into annotation_results_ch
    

    script:
    """
    mkdir annotation_results
    vep -i ${sample_id} --cache --dir_cache /home/huk/Documents/vep_cache/GRCh37/vep --offline --fasta /home/huk/Documents/vep_cache/GRCh37/Homo_sapiens.GRCh37.cdna.all.fa --hgvs --port 3337 -o ${sample_id.baseName}_vep.txt --everything --fork $task.cpus --plugin CADD --custom /home/huk/Documents/vep_cache/GRCh37/clinvar_20220528.vcf.gz,ClinVar,vcf,exact,0,CLNSIG
    """
}






