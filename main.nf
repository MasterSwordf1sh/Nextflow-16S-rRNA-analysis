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
Memory