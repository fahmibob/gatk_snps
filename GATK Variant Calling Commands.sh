#Alignment
bwa mem -R "@RG\tID:iriomote\tSM:061" -t 2 -K 100000000 ReferenceGenome.fasta repaired_trimmed_Bjaponica-061_S156_R1.fastq repaired_trimmed_Bjaponica-061_S156_R2.fastq -o aligned_061.sam

#MarkDuplicatesSpark
gatk MarkDuplicatesSpark -I aligned_061.sam -M 061_dedup_metrics.txt -O 061_sorted_dedup_reads.bam

#Collects
#CollectAlignmentSummaryMetrics 
java -jar picard.jar CollectAlignmentSummaryMetrics -R ReferenceGenome.fasta -I 061_sorted_dedup_reads.bam -O 061_alignment_metrics.txt
#CollectInsertSizeMetrics 
java -jar picard.jar CollectInsertSizeMetrics I=061_sorted_dedup_reads.bam O=061_insert_metrics.txt HISTOGRAM_FILE=061_insert_size_histogram.pdf
samtools depth -a 061_sorted_dedup_reads.bam > 061_depth_out.txt

##before call variant need reference genome index file
#Creating the FASTA sequence dictionary file
gatk CreateSequenceDictionary -R ReferenceGenome.fasta
#Creating the fasta index file
samtools faidx ReferenceGenome.fasta

#call_variant
gatk HaplotypeCaller -R ReferenceGenome.fasta -I 061_sorted_dedup_reads.bam -O 061_raw_variants.vcf

##extract_SNPs_indels
#extract SNPs
gatk SelectVariants -R ReferenceGenome.fasta -V 061_raw_variants.vcf -select-type SNP -O 061_raw_snps.vcf
#extract indels
gatk SelectVariants -R ReferenceGenome.fasta -V 061_raw_variants.vcf -select-type INDEL -O 061_raw_indels.vcf

##variantFiltration
#variantFiltration SNPs
gatk VariantFiltration \
        -R ReferenceGenome.fasta \
        -V 061_raw_snps.vcf \
        -O 061_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
#variantFiltration Indels
gatk VariantFiltration \
        -R ReferenceGenome.fasta \
        -V 061_raw_indels.vcf \
        -O 061_filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0"

##excludeFilteredVariants
#excludeFilteredVariants SNPs
gatk SelectVariants \
        --exclude-filtered \
        -V 061_filtered_snps.vcf \
        -O 061_bqsr_snps.vcf
#excludeFilteredVariants indels
gatk SelectVariants \
        --exclude-filtered \
        -V 061_raw_indels.vcf \
        -O 061_bqsr_indels.vcf

##Before BQSR
modify the READ GROUPS by AddOrReplaceReadGroups (Picard)
https://gatk.broadinstitute.org/hc/en-us/articles/360037226472
https://gatk.broadinstitute.org/hc/en-us/community/posts/360077416331-MGI-DNBSEQ-T7-DNA-Sequencer-Nebula-Is-it-Compatible-with-GATK-Best-Practice-Pipeline-
java -jar picard.jar AddOrReplaceReadGroups \
       I=061_sorted_dedup_reads.bam \
       O=ModifiedRG_061_sorted_dedup_reads.bam \
       RGID=iriomote \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=061
cek: samtools view -H ModifiedRG_061_sorted_dedup_reads.bam | grep '^@RG'

##bqsr_1 
#BaseRecalibrator
gatk BaseRecalibrator \
        -R ReferenceGenome.fasta \
        -I ModifiedRG_061_sorted_dedup_reads.bam \
        --known-sites 061_bqsr_snps.vcf \
        --known-sites 061_bqsr_indels.vcf \
        -O 061_recal_data.table
#ApplyBQSR
gatk ApplyBQSR \
        -R ReferenceGenome.fasta \
        -I ModifiedRG_061_sorted_dedup_reads.bam \
        -bqsr 061_recal_data.table \
        -O 061_recal_reads.bam

#bqsr_2
gatk BaseRecalibrator \
        -R ReferenceGenome.fasta \
        -I 061_recal_reads.bam \
        --known-sites 061_bqsr_snps.vcf \
        --known-sites 061_bqsr_indels.vcf \
        -O 061_post_recal_data.table

#call_variant_2
gatk HaplotypeCaller \
        -R ReferenceGenome.fasta \
        -I 061_recal_reads.bam \
        -O 061_raw_variants_recal.vcf

##extract_SNPs_indels2
#extract_SNPs2
gatk SelectVariants \
            -R ReferenceGenome.fasta \
            -V 061_raw_variants_recal.vcf \
            -select-type SNP \
            -O 061_raw_snps_recal.vcf
#extract_indels2
gatk SelectVariants \
            -R ReferenceGenome.fasta \
            -V 061_raw_variants_recal.vcf \
            -select-type INDEL \
            -O 061_raw_indels_recal.vcf

##variantFiltration2
#variantFiltration SNPs2
gatk VariantFiltration \
        -R ReferenceGenome.fasta \
        -V 061_raw_snps_recal.vcf \
        -O 061_filtered_snps_final.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
#variantFiltration Indels2
gatk VariantFiltration \
        -R ReferenceGenome.fasta \
        -V 061_raw_indels_recal.vcf \
        -O 061_filtered_indels_final.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0"

#compileStatistics
used bash parse_metrics.sh

