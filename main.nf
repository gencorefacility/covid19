/*  COVID-19 Pipeline
 *  Usage: nextflow run /path/to/main.nf
 *
 *  Author: Mohammed Khalfan < mkhalfan@nyu.edu >
 *  NYU Center for Genetics and System Biology 2020
 */

// Setting some defaults here,
// can be overridden in config or via command line
params.out = "${params.outdir}/out"
params.tmpdir = "${params.outdir}/gatk_temp"
params.snpEff_config = "static/snpEff.config"

// Define modules here
BWA = 'bwa/intel/0.7.17'
PICARD = 'picard/2.17.11'
GATK = 'gatk/4.1.3.0'
R = 'r/intel/3.4.2'
SAMTOOLS = 'samtools/intel/1.9'
TRIMMOMATIC = 'trimmomatic/0.36'
SNPEFF = 'snpeff/4.3i'
PYPAIRIX = 'pypairix/intel/0.2.3'
HTSLIB = 'htslib/intel/1.4.1'
DEEPTOOLS = 'deeptools/intel/2.4.2'
JVARKIT = 'jvarkit/base'
PYSAM = 'pysam/intel/python3.6/0.14.1'
PILON = 'pilon/1.23'

println "reads: $params.reads"
println "ref: $params.ref"
println "outdir: $params.out"

ref = file(params.ref)
snpeff_config = file(params.snpEff_config)
primers = file(params.primers)

// Prepare the fastq read pairs for input.
// Use the size parameter to not auto-group, and instead
// use the mapping through getBaseName() and subtract
// two regexs to get the ID.
// This enables support for CGSB sequence data file naming format
Channel
    .fromFilePairs( params.reads, size: -1)
    { file -> file.getBaseName() - ~/${params.read_pair_regex}/ - ~/.fastq/ }
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .set { read_pairs_ch }

process trim {
    publishDir "${params.out}/trimmed", mode:'copy'

    input:
    file(primers)
    set pair_id,
        file(reads) from read_pairs_ch

    output:
    set val(pair_id),
	file("${pair_id}_trimmed_1.fq.gz"),
	file("${pair_id}_trimmed_2.fq.gz") \
	into trimmed_ch

    script:
    """
    module load $TRIMMOMATIC
    java -jar \$TRIMMOMATIC_JAR \
	PE \
	-phred33 \
	-threads ${task.cpus} \
	${reads[0]} \
	${reads[1]} \
	${pair_id}_trimmed_1.fq.gz \
	${pair_id}.unpair_trimmed_1.fq.gz \
	${pair_id}_trimmed_2.fq.gz \
	${pair_id}.unpair_trimmed_2.fq.gz \
	ILLUMINACLIP:${params.adapters}:2:30:10:8:true \
        ILLUMINACLIP:${primers}:2:30:10:8:true \
	LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20 
    """
}

process align {
    publishDir "${params.out}/aligned_reads", mode:'copy'
	
    input:
    set pair_id, 
	file(read_1),
	file(read_2) from trimmed_ch
     
    output:
    val(pair_id) into jbrowse_pair_id_ch
    set val(sample_id), 
	file("${pair_id}_aligned_reads.bam") \
	into aligned_reads_ch
    set val(pair_id),
        file("${pair_id}_aligned_reads.bam"),
	file("${pair_id}_aligned_reads.bai") \
	into individual_bw_ch
	
    script:
    sample_id=pair_id - ~/${params.grouping_regex}/
    readGroup = "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${pair_id}"
    """
    module load $BWA
    bwa mem \
	-K 100000000 \
	-v 3 -t ${task.cpus} \
	-Y \
	-R \"${readGroup}\" \
	$ref \
	$read_1 \
	$read_2 \
	> ${pair_id}_aligned_reads.sam

    module load $PICARD
    java -jar \$PICARD_JAR SortSam \
	I=${pair_id}_aligned_reads.sam \
	O=${pair_id}_aligned_reads.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true
    """
}

process mergeBam{
    publishDir "${params.out}/merged", mode:'copy'
    input:
    set sample_id, file(sams) \
	from aligned_reads_ch
	.groupTuple(size: 2, remainder: true)

    output:
    set val(sample_id),
	file("{*fixed.bam,*unmerged.bam}") \
	into merged_bam_ch

    script:
    if( sams.size() == 2 )
    """
    module load $PICARD
    java -jar \$PICARD_JAR MergeSamFiles \
	I=${sams[0]} \
	I=${sams[1]} \
	O=${sample_id}_merged.bam
    java -jar \$PICARD_JAR AddOrReplaceReadGroups \
	I=${sample_id}_merged.bam \
        O=${sample_id}_fixed.bam \
        RGID=${sample_id} \
        RGLB=${sample_id} \
        RGPL=${params.pl} \
        RGPU=${params.fcid} \
        RGSM=${sample_id}
    """

    else
    """
    # need to 'output' the file so it goes in the channel
    mv ${sams[0]} ${sample_id}_unmerged.bam
    """

}

process markDuplicatesSpark  {
    publishDir "${params.out}/sorted", mode:'copy'

    input:
    set val(sample_id), 
	file(bam) from merged_bam_ch

    output:
    val(sample_id) into jbrowse_sample_id_ch
    set val(sample_id),
	file("${sample_id}_sorted_dedup.bam") \
	into sorted_dedup_bam_ch, sorted_dedup_ch_for_metrics, downsample_bam_ch, pilon_ch
    set val(sample_id),
        file("${sample_id}_sorted_dedup.bam"),
        file("${sample_id}_sorted_dedup.bam.bai") \
	into merged_bw_ch
    set val(sample_id),
	file("${sample_id}_dedup_metrics.txt") into dedup_qc_ch

    script:
    """
    module load $GATK
    mkdir -p $params.tmpdir/$workflow.runName/$sample_id
    gatk --java-options "-Djava.io.tmpdir=${params.tmpdir}/${workflow.runName}/${sample_id}" \
	 MarkDuplicatesSpark \
	-I ${bam} \
	-M ${sample_id}_dedup_metrics.txt \
	-O ${sample_id}_sorted_dedup.bam
    rm -r $params.tmpdir/$workflow.runName/$sample_id
    """ 
}

process getMetrics {
    publishDir "${params.out}/metrics", mode:'copy'

    input:
    set val(sample_id),
	file(sorted_dedup_reads) from sorted_dedup_ch_for_metrics

    output:
    set val(sample_id), 
            file("${sample_id}_alignment_metrics.txt"),
            file("${sample_id}_insert_metrics.txt"),
            file("${sample_id}_insert_size_histogram.pdf"),
            file("${sample_id}_depth_out.txt") \
            into metrics_output

    script:
    """
    module load $PICARD
    module load $R
    module load $SAMTOOLS
    java -jar \$PICARD_JAR \
        CollectAlignmentSummaryMetrics \
	R=${params.ref} \
        I=${sorted_dedup_reads} \
	O=${sample_id}_alignment_metrics.txt
    java -jar \$PICARD_JAR \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	OUTPUT=${sample_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${sample_id}_insert_size_histogram.pdf 
    samtools depth -a ${sorted_dedup_reads} > ${sample_id}_depth_out.txt
    """
}

process pilon{
    publishDir "${params.out}/pilon", mode:'copy'

    input:
    set val(sample_id),
	file(preprocessed_bam) from pilon_ch

    output:
    file("${sample_id}_pilon.vcf") into pilon_bzip_tabix_vcf_ch
    file '*' into pilon_out_ch

    script:
    """
    module load $PILON
    java -Xmx16G -jar \$PILON_JAR \
	--genome $ref \
	--bam $preprocessed_bam \
	--fix bases \
	--changes \
	--vcf \
	--threads ${task.cpus} \
	--mindepth 10 \
	--output ${sample_id}_pilon_g
    
    module load $GATK
    gatk SelectVariants \
	-V ${sample_id}_pilon_g.vcf \
	-O ${sample_id}_pilon.vcf \
	--exclude-non-variants \
	--exclude-filtered
    """
}

process haplotypeCaller {
    input:
    set val(sample_id), 
	file(preprocessed_bam) from sorted_dedup_bam_ch

    output:
    set val(sample_id), 
	file("${sample_id}_raw_variants.vcf") into hc_output_ch

    script:
    """
    module load $GATK
    gatk HaplotypeCaller \
	-R $ref \
	-I $preprocessed_bam \
	-O ${sample_id}_raw_variants.vcf \
	-ploidy 1	
    """
}

process selectVariants {
    input:
    set val(sample_id), 
	file(raw_variants) from hc_output_ch

    output:
    set val(sample_id),
	file("${sample_id}_raw_snps.vcf") \
	into raw_snps_ch, raw_snps_qc_ch
    set val(sample_id),
	file("${sample_id}_raw_indels.vcf") into raw_indels_ch

    script:
    """
    module load $GATK
    gatk SelectVariants \
	-R $ref \
	-V $raw_variants \
	-select-type SNP \
	-O ${sample_id}_raw_snps.vcf
    gatk SelectVariants \
        -R $ref \
        -V $raw_variants \
        -select-type INDEL \
        -O ${sample_id}_raw_indels.vcf
    """
}

process filterSnps {
    publishDir "${params.out}/filtered_snps", mode:'copy'
    
    input:
    set val(sample_id), 
	file(raw_snps) from raw_snps_ch

    output:
    set val(sample_id),
        file("${sample_id}_filtered_snps.vcf") \
        into filtered_snps_qc_ch
    set val(sample_id),
	file("${sample_id}_filtered_snps_eaf.vcf") \
	into snpeff_ch
    set val(sample_id),
        file("${sample_id}_consensus_snps.vcf") \
        into consensus_snps_ch
    file("${sample_id}_consensus_snps.vcf") \
	into cons_bzip_tabix_vcf_ch

    script:
    """
    module load $GATK
    gatk VariantFiltration \
	-R $ref \
	-V $raw_snps \
	-O ${sample_id}_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

    # This script generates the _consensus_snps.vcf
    # and _eaf.vcf (empirical AF) files
    module load $PYSAM
    filter_variants.py ${sample_id}
    """
}

process filterIndels {
    publishDir "${params.out}/filtered_indels", mode:'copy'
    input:
    set val(sample_id),
	file(raw_indels) from raw_indels_ch

    output:
    file("${sample_id}_filtered_indels.vcf") into indel_bzip_tabix_vcf_ch

    script:
    """
    module load $GATK
    gatk VariantFiltration \
        -R $ref \
        -V $raw_indels \
        -O ${sample_id}_filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}

process consensus {
    publishDir "${params.out}/consensus", mode:'copy' 

    input:
    set val(sample_id), 
	file(filtered_snps) \
	from consensus_snps_ch

    output:
    file("${sample_id}.fasta") into consensus_ch

    script:
    """
    module load $GATK
    gatk IndexFeatureFile \
	-F $filtered_snps
    gatk FastaAlternateReferenceMaker \
	-R $ref \
	-O ${sample_id}.fasta \
	-V $filtered_snps
    """
}

process snpEff{
    publishDir "${params.out}/snpEff", mode:'copy'

    input:
    file(snpeff_config)
    set val(sample_id), 
	file(snps) \
	from snpeff_ch

    output:
    file '*' into snpeff_out
    file("${sample_id}_filtered_snps.ann.vcf") into snpeff_bzip_tabix_vcf_ch

    script:
    """
    module load $SNPEFF
    java -jar \$SNPEFF_JAR -v \
        -c $snpeff_config \
        SARS-CoV2_NC_045512.2 \
        $snps > ${sample_id}_filtered_snps.ann.vcf
    """
}

process make_bw{
    publishDir "${params.out}/bigwig", mode:'copy'

    input:
    /* id can be sample_id or pair_id */
    set val(id), 
	file(bam),
	file(bam_index) \
	from individual_bw_ch.mix(merged_bw_ch)

    output:
    file("${id}_coverage.bam.bw") into jbrowse_bw_ch 

    script:
    """
    module load $DEEPTOOLS
    bamCoverage \
        -p max  \
        --bam $bam \
        -o ${id}_coverage.bam.bw
    """
}

process downsample_bam{
    input:
    set val(sample_id), file(bam) from downsample_bam_ch

    output:
    set file("${sample_id}_downsampled.bam"),
        file("${sample_id}_downsampled.bam.bai") into jbrowse_bam_ch

    script:
    """
    module load $JVARKIT
    module load $SAMTOOLS
    java -jar \$SORTSAMREFNAME_JAR \
        --samoutputformat BAM \
        $bam |\
        java -jar \$BIOSTAR_JAR \
        -n 75 \
        --samoutputformat BAM |\
        samtools sort -o ${sample_id}_downsampled.bam
    samtools index ${sample_id}_downsampled.bam
    """
}

process bzip_tabix_vcf{
    input:
    file(vcf) from pilon_bzip_tabix_vcf_ch
	.mix(cons_bzip_tabix_vcf_ch)
	.mix(indel_bzip_tabix_vcf_ch)
	.mix(snpeff_bzip_tabix_vcf_ch)

    output:
    file("*.vcf.gz*") into jbrowse_vcf_ch

    script:
    """
    module load $HTSLIB
    module load $PYPAIRIX
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    """
}

process jbrowse{
    publishDir "${params.out}/trackList", mode:'copy'

    input:
    val pair_ids from jbrowse_pair_id_ch.collect()
    val sample_ids from jbrowse_sample_id_ch.collect()
    file '*' from jbrowse_bw_ch.collect()
    file '*' from jbrowse_bam_ch.collect()
    file '*' from jbrowse_vcf_ch.collect()

    output:
    file '*.json' into trackList_ch

    when:
    false

    script:
    """
    copy-data-to-jbrowse.sh "${pair_ids}" "${sample_ids}" $params.fcid "$params.grouping_regex"
    """
}

process qc {
    input:
    set val(sample_id),
	file("${sample_id}_alignment_metrics.txt"),
	file("${sample_id}_insert_metrics.txt"),
	file("${sample_id}_insert_size_histogram.pdf"),
	file("${sample_id}_depth_out.txt"),
	file("${sample_id}_dedup_metrics.txt"),
	file("${sample_id}_raw_snps.vcf"),
        file("${sample_id}_filtered_snps.vcf") \
	from metrics_output
	.join(dedup_qc_ch)
	.join(raw_snps_qc_ch)
	.join(filtered_snps_qc_ch)

    output:
    file("${sample_id}_report.csv") into parse_metrics_output

    script:
    """
    parse_metrics.sh ${sample_id} > ${sample_id}_report.csv 
    """
}

/* Process qc above creates a report for each sample.
 * Below we compile these into a single report.
 */
parse_metrics_output.collectFile(name: "${workflow.runName}_report.csv", keepHeader: true, storeDir: "${params.out}/reports")

