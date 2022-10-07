nextflow.enable.dsl=2
params.inputsam="aligned_061.sam"
params.ref="ReferenceGenome.fasta"
params.index="$baseDir/index/*"
params.forking=4
params.job="default"
params.parse="$baseDir/bin/parse_metrics.sh"

workflow {
  if (params.inputsam) {
    Channel.fromPath(params.inputsam)
      .map {file -> tuple(file.simpleName, file)}
      .set{pairs_CH}
    //Channel.fromPath(params.index)
    //  .collect()
    //  .set{index_CH}

    Channel.fromPath(params.ref)
      .set{ref_CH}
    Channel.fromPath(params.index)
      .collect()
      .set{index_CH}

    //bwa_mem(pairs_CH, index_CH, ref_CH)
    markdupplicate(pairs_CH)
    collectAlignment(markdupplicate.out, ref_CH)
    call_variant(markdupplicate.out, ref_CH, index_CH)
    extract_SNPs_indels(call_variant.out, ref_CH, index_CH)
    variantFiltration(extract_SNPs_indels.out, ref_CH, index_CH)
    excludeFilteredVariants(variantFiltration.out, ref_CH, index_CH)
    bqsr_1(excludeFilteredVariants.out, markdupplicate.out, ref_CH, index_CH)
    bqsr_2(excludeFilteredVariants.out, bqsr_1.out, ref_CH, index_CH)
    //analyzeCovariates(bqsr_1.out, bqsr_2.out)
    call_variant_2(bqsr_1.out, ref_CH, index_CH)
    extract_SNPs_indels2(call_variant_2.out, ref_CH, index_CH)
    variantFiltration2(extract_SNPs_indels2.out, ref_CH, index_CH)
    compileStatistics(collectAlignment.out, markdupplicate.out, extract_SNPs_indels.out, variantFiltration.out, extract_SNPs_indels2.out, variantFiltration2.out)
  }

}

//process bwa_mem {
  //publishDir "$baseDir/output/$params.job/1_bwamem", mode: 'copy'
  //maxForks params.forking
  //errorStrategy 'ignore'

  //input:
  //tuple val(sampleName), path(forward), path(reverse)
  //path(index)
  //path(reference)

  //output:
  //tuple val(sampleName), path('*.sam')

  //script:
  //"""

  //header=$(zcat $1 | head -n 1)
  //id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
  //sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
  //echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"
  //\$(echo "@RG\tID:\$id\tSM:\$id"_"\$sm\tLB:\$id"_"\$sm\tPL:ILLUMINA")
  //  bwa mem -t 2\
  //  -K 100000000 \
  //  $reference \
  //  $forward \
  //  $reverse \
  //  > aligned_$sampleName'.sam'
  //"""
//}

process markdupplicate {
  publishDir "$baseDir/output/$params.job/2_markdupplicate", mode: 'copy'
  maxForks params.forking


  input:
  tuple val(sampleName), path(samfile)

  output:
  tuple val(sampleName), path('*dedup_metrics.txt'), path('*sorted_dedup_reads.bam')

  script:
  """
    java -jar ~/picard.jar SortSam I=$samfile O=sorted$samfile SORT_ORDER=coordinate
    gatk MarkDuplicates \
        -I sorted$samfile \
        -M $sampleName'dedup_metrics.txt' \
        -O $sampleName'sorted_dedup_reads.bam'
  """
}

process collectAlignment {
  publishDir "$baseDir/output/$params.job/3_collectAlignment", mode: 'copy'
  maxForks params.forking
  errorStrategy 'ignore'

  input:
  tuple val(sampleName), path(dedup), path(sorted_dedup)
  path(reference)

  output:
  tuple path('*alignment_metrics.txt'), path('*insert_metrics.txt'), path('*insert_size_histogram.pdf'), path('*depth_out.txt')

  script:
  """
    java -jar ~/picard.jar \
        CollectAlignmentSummaryMetrics \
        -R $reference \
        -I $sorted_dedup \
        -O $sampleName'alignment_metrics.txt'\

    java -jar ~/picard.jar \
        CollectInsertSizeMetrics \
        I=$sorted_dedup \
        O=$sampleName'insert_metrics.txt' \
        HISTOGRAM_FILE=$sampleName'insert_size_histogram.pdf'

    samtools depth -a $sorted_dedup > $sampleName'depth_out.txt'
  """
}

process call_variant {
  publishDir "$baseDir/output/$params.job/4_call_variant", mode: 'copy'
  maxForks params.forking
  errorStrategy 'ignore'

  input:
  tuple val(sampleName), path(dedup), path(sorted_dedup)
  path(reference)
  path(index)

  output:
  tuple val(sampleName), path('*raw_variants.vcf')


  script:
  """
    samtools index $sorted_dedup
    gatk HaplotypeCaller \
        -R $reference \
        -I $sorted_dedup \
        -O $sampleName'raw_variants.vcf'
  """
}

process extract_SNPs_indels {
    publishDir "$baseDir/output/$params.job/5_extract_SNPs_indels", mode: 'copy'
    maxForks params.forking
    errorStrategy 'ignore'

    input:
    tuple val(sampleName), path(rawvariants)
    path(reference)
    path(index)

    output:
    tuple val(sampleName), path('*raw_snps.vcf'), path('*raw_indels.vcf')

    script:
    """
        gatk SelectVariants \
            -R $reference \
            -V $rawvariants \
            -select-type SNP \
            -O $sampleName'raw_snps.vcf'

        gatk SelectVariants \
            -R $reference \
            -V $rawvariants \
            -select-type INDEL \
            -O $sampleName'raw_indels.vcf'
    """
}

process variantFiltration  {
    publishDir "$baseDir/output/$params.job/6_variantFiltration", mode: 'copy'
    maxForks params.forking
    errorStrategy 'ignore'

    input:
    tuple val(sampleName), path(raw_snps), path(raw_indels)
    path(reference)
    path(index)

    output:
    tuple val(sampleName), path('*filtered_snps.vcf'), path('*filtered_indels.vcf')

    script:
    """
        gatk VariantFiltration \
        -R $reference \
        -V $raw_snps \
        -O $sampleName'filtered_snps.vcf' \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

         gatk VariantFiltration \
        -R $reference \
        -V $raw_indels \
        -O $sampleName'filtered_indels.vcf' \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}


process excludeFilteredVariants  {
    publishDir "$baseDir/output/$params.job/7_excludeFilteredVariants", mode: 'copy'
    maxForks params.forking
    errorStrategy 'ignore'

    input:
    tuple val(sampleName), path(filtered_snps), path(filtered_indels)
    path(reference)
    path(index)

    output:
    tuple val(sampleName), path('*bqsr_snps.vcf'), path('*bqsr_indels.vcf')

    script:
    """
        gatk SelectVariants \
        --exclude-filtered \
        -V $filtered_snps \
        -O $sampleName'bqsr_snps.vcf'

        gatk SelectVariants \
        --exclude-filtered \
        -V $filtered_indels \
        -O $sampleName'bqsr_indels.vcf'
    """
}

process bqsr_1  {
    publishDir "$baseDir/output/$params.job/8_bqsr_1", mode: 'copy'
    maxForks params.forking
    errorStrategy 'ignore'

    input:
    tuple val(sampleName), path(bqsr_snps), path(bqsr_indels)
    tuple val(sampleName2), path(dedup), path(sorted_dedup)
    path(reference)
    path(index)

    output:
    tuple val(sampleName), path('*recal_data.table'), path('*recal_reads.bam')

    script:
    """
    java -jar ~/picard.jar AddOrReplaceReadGroups \
      I=$sorted_dedup \
      O=ModifiedRG$sorted_dedup \
      RGID=iriomote \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=unit1 \
      RGSM=061

    gatk IndexFeatureFile -I $bqsr_snps
    gatk IndexFeatureFile -I $bqsr_indels


        gatk BaseRecalibrator \
        -R $reference \
        -I ModifiedRG$sorted_dedup \
        --known-sites $bqsr_snps \
        --known-sites $bqsr_indels \
        -O $sampleName'recal_data.table'

        gatk ApplyBQSR \
        -R $reference \
        -I ModifiedRG$sorted_dedup \
        -bqsr $sampleName'recal_data.table' \
        -O $sampleName'recal_reads.bam' \
    """
}

process bqsr_2  {
    publishDir "$baseDir/output/$params.job/9_bqsr_2", mode: 'copy'
    maxForks params.forking
    errorStrategy 'ignore'

    input:
    tuple val(sampleName), path(bqsr_snps), path(bqsr_indels)
    tuple val(sampleName2), path(recal_data), path(recal_reads)
    path(reference)
    path(index)

    output:
    path('*post_recal_data.table')

    script:
    """
    gatk IndexFeatureFile -I $bqsr_snps
    gatk IndexFeatureFile -I $bqsr_indels
        gatk BaseRecalibrator \
        -R $reference \
        -I $recal_reads \
        --known-sites $bqsr_snps \
        --known-sites $bqsr_indels \
        -O $sampleName'post_recal_data.table'
    """
}

//process analyzeCovariates {
//    publishDir "$baseDir/output/$params.job/10_analyzeCovariates", mode: 'copy'
//    maxForks params.forking
//    errorStrategy 'ignore'

//    input:
//    tuple val(sampleName), path(recal_data), path(recal_reads)
//    path(post_recal_data)

//    output:
//    tuple val(sampleName), path('*recalibration_plots.pdf')

//    script:
  //  """
    //  gatk AnalyzeCovariates \
    //  -before $recal_data \
    //  -after $post_recal_data \
    //  -plots $sampleName'recalibration_plots.pdf'
    //"""
//}

process call_variant_2 {
  publishDir "$baseDir/output/$params.job/11_call_variant_2", mode: 'copy'
  maxForks params.forking
  errorStrategy 'ignore'

  input:
  tuple val(sampleName), path(recal_data), path(recal_reads)
  path(reference)
  path(index)

  output:
  tuple val(sampleName), path('*raw_variants_recal.vcf')

  script:
  """
    samtools index $recal_reads
    gatk HaplotypeCaller \
        -R $reference \
        -I $recal_reads \
        -O $sampleName'raw_variants_recal.vcf'
  """
}

process extract_SNPs_indels2 {
    publishDir "$baseDir/output/$params.job/12_extract_SNPs_indels2", mode: 'copy'
    maxForks params.forking
    errorStrategy 'ignore'

    input:
    tuple val(sampleName), path(rawvariantsrecal)
    path(reference)
    path(index)

    output:
    tuple val(sampleName), path('*raw_snps_recal.vcf'), path('*raw_indels_recal.vcf')

    script:
    """
        gatk SelectVariants \
            -R $reference \
            -V $rawvariantsrecal \
            -select-type SNP \
            -O $sampleName'raw_snps_recal.vcf'

        gatk SelectVariants \
            -R $reference \
            -V $rawvariantsrecal \
            -select-type INDEL \
            -O $sampleName'raw_indels_recal.vcf'
    """
}

process variantFiltration2  {
    publishDir "$baseDir/output/$params.job/13_variantFiltration2", mode: 'copy'
    maxForks params.forking
    errorStrategy 'ignore'

    input:
    tuple val(sampleName), path(raw_snps_recal), path(raw_indels_recal)
    path(reference)
    path(index)

    output:
    tuple val(sampleName), path('*filtered_snps_final.vcf'), path('*filtered_indels_final.vcf')

    script:
    """
        gatk VariantFiltration \
        -R $reference \
        -V $raw_snps_recal \
        -O $sampleName'filtered_snps_final.vcf' \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

         gatk VariantFiltration \
        -R $reference \
        -V $raw_indels_recal \
        -O $sampleName'filtered_indels_final.vcf' \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}

process compileStatistics {
    publishDir "$baseDir/output/$params.job/14_compileStatistics", mode: 'copy'
    maxForks params.forking
    errorStrategy 'ignore'

    input:
    //collectAlignment
    tuple path(alignment_metrics), path(insert_metrics), path(insert_size_histogram), path(depth_out)
    //markdupplicate
    tuple val(sampleName), path(dedup_metrics), path(sorted_dedup_reads)
    //extract_SNPs
    tuple val(sampleName2), path(raw_snps), path(raw_indels)
    //variantFiltration
    tuple val(sampleName3), path(filtered_snps), path(filtered_indels)
    //extract_SNPs_2
    tuple val(sampleName4), path(raw_snps_recal), path(raw_indels_recal)
    //variantFiltration
    tuple val(sampleName5), path(filtered_snps_final), path(filtered_indels_final)


    output:
    tuple val(sampleName), path('*filtered_snps_final.vcf'), path('*filtered_indels_final.vcf')

    script:
    parse="$baseDir/bin/parse_metrics.sh"
    """
        bash $params.parse $alignment_metrics \
            $insert_metrics \
            $dedup_metrics \
            $filtered_snps\
            $filtered_snps_final \
            $depth_out \
            $sampleName > $sampleName'_report.csv'
    """
}
