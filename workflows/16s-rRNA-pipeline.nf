#!/usr/bin/env nextflow

/*
 * 16S rRNA Analysis Pipeline
 * A complete workflow for processing and analyzing 16S rRNA amplicon sequencing data
 */

// Pipeline parameter defaults
params {
    // Input
    reads = "data/raw/*_R{1,2}.fastq.gz"
    metadata = "data/metadata.tsv"
    
    // Output directories
    outdir = "results"
    
    // Resources
    cpus = 8
    memory = '16.GB'
    
    // Databases
    silva_db = "databases/silva_138_99_16S.qza"
    
    // Trimming parameters
    trim_front_f = 0
    trim_front_r = 0
    trunc_len_f = 0      // 0 means no truncation
    trunc_len_r = 0      // 0 means no truncation
    min_len = 100
    max_ee = 2.0
    
    // ASV inference parameters
    pooling = "pseudo"   // Options: pseudo, independent, or none
    
    // Taxonomy classification parameters
    confidence = 0.8
    
    // Alpha diversity metrics
    alpha_metrics = "shannon,observed_features,faith_pd,evenness"
    
    // Beta diversity metrics
    beta_metrics = "jaccard,bray_curtis,unweighted_unifrac,weighted_unifrac"
    
    // Visualization
    n_taxa_barplot = 10  // Number of top taxa to show in barplots
    
    // Skip options
    skip_fastqc = false
    skip_multiqc = false
}

// Print pipeline info
log.info """
=============================================================
16S rRNA ANALYSIS PIPELINE
=============================================================
Input reads       : ${params.reads}
Metadata file     : ${params.metadata}
Output directory  : ${params.outdir}
SILVA database    : ${params.silva_db}
CPUs per task     : ${params.cpus}
Memory per task   : ${params.memory}
=============================================================
"""

// Define channels from input files
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { raw_reads_ch }

// Step 1: Quality control with FastQC
process FASTQC {
    tag "FastQC on ${sample_id}"
    publishDir "${params.outdir}/1_fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from raw_reads_ch
    
    output:
    path("fastqc_${sample_id}_logs") into fastqc_ch
    tuple val(sample_id), path(reads) into trim_ch
    
    when:
    !params.skip_fastqc
    
    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -q -t ${params.cpus} ${reads[0]} ${reads[1]}
    """
}

// Step 2: Adapter trimming and quality filtering with Trimmomatic
process TRIMMOMATIC {
    tag "Trimmomatic on ${sample_id}"
    publishDir "${params.outdir}/2_trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from trim_ch
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz") into trimmed_reads_ch
    
    script:
    """
    trimmomatic PE -threads ${params.cpus} \
        ${reads[0]} ${reads[1]} \
        ${sample_id}_trimmed_R1.fastq.gz ${sample_id}_unpaired_R1.fastq.gz \
        ${sample_id}_trimmed_R2.fastq.gz ${sample_id}_unpaired_R2.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:${params.min_len}
    """
}

// Step 3: Import data to QIIME2 format
process QIIME2_IMPORT {
    tag "QIIME2 import on ${sample_id}"
    publishDir "${params.outdir}/3_qiime2_import", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from trimmed_reads_ch
    
    output:
    tuple val(sample_id), path("${sample_id}.qza") into demux_ch
    
    script:
    """
    # Create manifest file
    echo "sample-id,absolute-filepath,direction" > manifest.csv
    echo "${sample_id},\$PWD/${reads[0]},forward" >> manifest.csv
    echo "${sample_id},\$PWD/${reads[1]},reverse" >> manifest.csv
    
    # Import data
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path manifest.csv \
        --output-path ${sample_id}.qza \
        --input-format PairedEndFastqManifestPhred33
    """
}

// Step 4: DADA2 denoising
process DADA2_DENOISE {
    tag "DADA2 denoising"
    publishDir "${params.outdir}/4_dada2", mode: 'copy'
    memory { 2.GB * task.cpus }
    
    input:
    tuple val(sample_id), path(demux) from demux_ch.collect()
    
    output:
    path("table.qza") into table_ch
    path("rep-seqs.qza") into repseqs_ch
    path("denoising-stats.qza") into stats_ch
    
    script:
    // Join all sample IDs with commas for the --i-demultiplexed-seqs parameter
    def demux_inputs = demux.collect { "--i-demultiplexed-seqs ${it}" }.join(' ')
    
    """
    # Merge all demux QZA files (if multiple samples)
    qiime feature-table merge-seqs \
        ${demux_inputs} \
        --o-merged-data merged-demux.qza
    
    # Run DADA2
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs merged-demux.qza \
        --p-trim-left-f ${params.trim_front_f} \
        --p-trim-left-r ${params.trim_front_r} \
        --p-trunc-len-f ${params.trunc_len_f} \
        --p-trunc-len-r ${params.trunc_len_r} \
        --p-max-ee ${params.max_ee} \
        --p-n-threads ${params.cpus} \
        --p-pooling-method ${params.pooling} \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza
    """
}

// Step 5: Generate feature table summaries
process FEATURE_TABLE_SUMMARY {
    tag "Feature table summary"
    publishDir "${params.outdir}/5_summaries", mode: 'copy'
    
    input:
    path(table) from table_ch
    path(repseqs) from repseqs_ch
    
    output:
    path("table-summary.qzv") into table_summary_ch
    path("rep-seqs-summary.qzv") into repseqs_summary_ch
    
    script:
    """
    qiime feature-table summarize \
        --i-table ${table} \
        --o-visualization table-summary.qzv
    
    qiime feature-table tabulate-seqs \
        --i-data ${repseqs} \
        --o-visualization rep-seqs-summary.qzv
    """
}

// Step 6: Phylogenetic tree construction
process PHYLOGENETIC_TREE {
    tag "Building phylogenetic tree"
    publishDir "${params.outdir}/6_phylogeny", mode: 'copy'
    
    input:
    path(repseqs) from repseqs_ch
    
    output:
    path("aligned-rep-seqs.qza") into aligned_seqs_ch
    path("masked-aligned-rep-seqs.qza") into masked_seqs_ch
    path("unrooted-tree.qza") into unrooted_tree_ch
    path("rooted-tree.qza") into rooted_tree_ch
    
    script:
    """
    # Sequence alignment
    qiime alignment mafft \
        --i-sequences ${repseqs} \
        --o-alignment aligned-rep-seqs.qza
    
    # Mask highly variable positions
    qiime alignment mask \
        --i-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza
    
    # Build unrooted tree
    qiime phylogeny fasttree \
        --i-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza
    
    # Root the tree at midpoint
    qiime phylogeny midpoint-root \
        --i-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza
    """
}

// Step 7: Taxonomic classification
process TAXONOMIC_CLASSIFICATION {
    tag "Taxonomic classification"
    publishDir "${params.outdir}/7_taxonomy", mode: 'copy'
    
    input:
    path(repseqs) from repseqs_ch
    
    output:
    path("taxonomy.qza") into taxonomy_ch
    path("taxonomy.qzv") into taxonomy_viz_ch
    
    script:
    """
    qiime feature-classifier classify-sklearn \
        --i-classifier ${params.silva_db} \
        --i-reads ${repseqs} \
        --p-confidence ${params.confidence} \
        --p-n-jobs ${params.cpus} \
        --o-classification taxonomy.qza
    
    qiime metadata tabulate \
        --m-input-file taxonomy.qza \
        --o-visualization taxonomy.qzv
    """
}

// Step 8: Taxonomic barplots
process TAXONOMIC_BARPLOTS {
    tag "Generating taxonomic barplots"
    publishDir "${params.outdir}/8_barplots", mode: 'copy'
    
    input:
    path(table) from table_ch
    path(taxonomy) from taxonomy_ch
    
    output:
    path("taxa-bar-plots.qzv") into barplots_ch
    
    script:
    """
    qiime taxa barplot \
        --i-table ${table} \
        --i-taxonomy ${taxonomy} \
        --m-metadata-file ${params.metadata} \
        --o-visualization taxa-bar-plots.qzv
    """
}

// Step 9: Alpha rarefaction curves
process ALPHA_RAREFACTION {
    tag "Alpha rarefaction curves"
    publishDir "${params.outdir}/9_diversity", mode: 'copy'
    
    input:
    path(table) from table_ch
    path(rooted_tree) from rooted_tree_ch
    
    output:
    path("alpha-rarefaction.qzv") into alpha_rarefaction_ch
    
    script:
    """
    qiime diversity alpha-rarefaction \
        --i-table ${table} \
        --i-phylogeny ${rooted_tree} \
        --p-max-depth auto \
        --m-metadata-file ${params.metadata} \
        --o-visualization alpha-rarefaction.qzv
    """
}

// Step 10: Core diversity metrics
process CORE_DIVERSITY_METRICS {
    tag "Core diversity metrics"
    publishDir "${params.outdir}/9_diversity", mode: 'copy'
    
    input:
    path(table) from table_ch
    path(rooted_tree) from rooted_tree_ch
    
    output:
    path("core-metrics-results/*") into diversity_results_ch
    
    script:
    """
    # Determine sampling depth from feature table
    SAMPLING_DEPTH=\$(qiime tools export \
        --input-path ${table} \
        --output-path table_export && \
        biom summarize-table \
        -i table_export/feature-table.biom | \
        grep -A 10 "Counts/sample detail" | \
        head -n 11 | \
        tail -n 1 | \
        awk '{print \$2}' | \
        sed 's/\\.0//')
    
    # Use 90% of the minimum sampling depth
    SAMPLING_DEPTH=\$(echo "\$SAMPLING_DEPTH * 0.9" | bc | sed 's/\\..*//g')
    
    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny ${rooted_tree} \
        --i-table ${table} \
        --p-sampling-depth \$SAMPLING_DEPTH \
        --m-metadata-file ${params.metadata} \
        --o-rarefied-table core-metrics-results/rarefied_table.qza \
        --o-faith-pd-vector core-metrics-results/faith_pd_vector.qza \
        --o-observed-features-vector core-metrics-results/observed_features_vector.qza \
        --o-shannon-vector core-metrics-results/shannon_vector.qza \
        --o-evenness-vector core-metrics-results/evenness_vector.qza \
        --o-unweighted-unifrac-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
        --o-weighted-unifrac-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
        --o-jaccard-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
        --o-bray-curtis-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
        --o-unweighted-unifrac-pcoa-results core-metrics-results/unweighted_unifrac_pcoa_results.qza \
        --o-weighted-unifrac-pcoa-results core-metrics-results/weighted_unifrac_pcoa_results.qza \
        --o-jaccard-pcoa-results core-metrics-results/jaccard_pcoa_results.qza \
        --o-bray-curtis-pcoa-results core-metrics-results/bray_curtis_pcoa_results.qza \
        --o-unweighted-unifrac-emperor core-metrics-results/unweighted_unifrac_emperor.qzv \
        --o-weighted-unifrac-emperor core-metrics-results/weighted_unifrac_emperor.qzv \
        --o-jaccard-emperor core-metrics-results/jaccard_emperor.qzv \
        --o-bray-curtis-emperor core-metrics-results/bray_curtis_emperor.qzv
    """
}

// Step 11: Alpha diversity statistical tests
process ALPHA_DIVERSITY_STATS {
    tag "Alpha diversity statistical analysis"
    publishDir "${params.outdir}/10_stats", mode: 'copy'
    
    input:
    path('core-metrics-results/*') from diversity_results_ch
    
    output:
    path("alpha-diversity-stats/*") into alpha_stats_ch
    
    script:
    def alpha_metrics = params.alpha_metrics.split(',')
    def alpha_cmds = alpha_metrics.collect { metric ->
        """
        qiime diversity alpha-group-significance \
            --i-alpha-diversity core-metrics-results/${metric}_vector.qza \
            --m-metadata-file ${params.metadata} \
            --o-visualization alpha-diversity-stats/${metric}-group-significance.qzv
        """
    }.join('\n')
    
    """
    mkdir -p alpha-diversity-stats
    
    ${alpha_cmds}
    """
}

// Step 12: Beta diversity statistical tests
process BETA_DIVERSITY_STATS {
    tag "Beta diversity statistical analysis"
    publishDir "${params.outdir}/10_stats", mode: 'copy'
    
    input:
    path('core-metrics-results/*') from diversity_results_ch
    
    output:
    path("beta-diversity-stats/*") into beta_stats_ch
    
    script:
    def beta_metrics = params.beta_metrics.split(',')
    def beta_cmds = beta_metrics.collect { metric ->
        """
        qiime diversity beta-group-significance \
            --i-distance-matrix core-metrics-results/${metric}_distance_matrix.qza \
            --m-metadata-file ${params.metadata} \
            --m-metadata-column sample-type \
            --p-permutations 999 \
            --o-visualization beta-diversity-stats/${metric}-group-significance.qzv
        
        qiime diversity adonis \
            --i-distance-matrix core-metrics-results/${metric}_distance_matrix.qza \
            --m-metadata-file ${params.metadata} \
            --p-formula "sample-type+treatment" \
            --p-permutations 999 \
            --o-visualization beta-diversity-stats/${metric}-adonis.qzv
        """
    }.join('\n')
    
    """
    mkdir -p beta-diversity-stats
    
    ${beta_cmds}
    """
}

// Step 13: ANCOM Differential Abundance Testing
process ANCOM_DIFFERENTIAL_ABUNDANCE {
    tag "ANCOM differential abundance"
    publishDir "${params.outdir}/11_differential_abundance", mode: 'copy'
    
    input:
    path(table) from table_ch
    path(taxonomy) from taxonomy_ch
    
    output:
    path("ancom-results/*") into ancom_results_ch
    
    script:
    """
    mkdir -p ancom-results
    
    # Collapse table to genus level
    qiime taxa collapse \
        --i-table ${table} \
        --i-taxonomy ${taxonomy} \
        --p-level 6 \
        --o-collapsed-table ancom-results/genus-table.qza
    
    # Filter low-frequency features
    qiime feature-table filter-features \
        --i-table ancom-results/genus-table.qza \
        --p-min-frequency 10 \
        --p-min-samples 2 \
        --o-filtered-table ancom-results/genus-table-filtered.qza
    
    # Add pseudocount for ANCOM
    qiime composition add-pseudocount \
        --i-table ancom-results/genus-table-filtered.qza \
        --o-composition-table ancom-results/genus-comp-table.qza
    
    # Run ANCOM
    qiime composition ancom \
        --i-table ancom-results/genus-comp-table.qza \
        --m-metadata-file ${params.metadata} \
        --m-metadata-column sample-type \
        --o-visualization ancom-results/ancom-sample-type.qzv
    
    # Also try with different columns if they exist
    if grep -q "treatment" ${params.metadata}; then
        qiime composition ancom \
            --i-table ancom-results/genus-comp-table.qza \
            --m-metadata-file ${params.metadata} \
            --m-metadata-column treatment \
            --o-visualization ancom-results/ancom-treatment.qzv
    fi
    """
}

// Step 14: Prepare report with MultiQC
process MULTIQC {
    tag "MultiQC report"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path('fastqc/*') from fastqc_ch.collect().ifEmpty([])
    path('stats/*') from stats_ch.ifEmpty([])
    
    output:
    path("multiqc_report.html") into multiqc_report_ch
    
    when:
    !params.skip_multiqc
    
    script:
    """
    multiqc .
    """
}

// Workflow completion message
workflow.onComplete {
    log.info """
    =============================================================
    16S rRNA ANALYSIS PIPELINE COMPLETED
    =============================================================
    Pipeline execution summary:
    ---------------------------
    Completed at  : ${workflow.complete}
    Duration      : ${workflow.duration}
    Success       : ${workflow.success}
    Work directory: ${workflow.workDir}
    Exit status   : ${workflow.exitStatus}
    =============================================================
    Results are available in: ${params.outdir}
    =============================================================
    """
}
