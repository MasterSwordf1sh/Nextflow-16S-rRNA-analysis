# 16S rRNA Analysis Pipeline

A complete Nextflow workflow for processing and analyzing 16S rRNA amplicon sequencing data from raw reads to statistical analysis.

![Pipeline Overview](docs/images/Nextflow_Pipeline.png)

## Overview

This pipeline automates the entire 16S rRNA analysis workflow, using QIIME2 and related tools to process raw sequencing data and generate comprehensive analysis results including taxonomic classification, diversity metrics, and differential abundance testing.

## Features

- Quality control and preprocessing of raw reads
- ASV inference with DADA2
- Taxonomic classification using the SILVA database
- Phylogenetic tree construction
- Alpha and beta diversity analysis
- Statistical testing for differential abundance (ANCOM)
- Comprehensive visualization outputs

## Requirements

- [Nextflow](https://www.nextflow.io/) (v21.04.0 or later)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) (recommended for reproducibility)
- [QIIME2](https://qiime2.org/) (2022.2 or later)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [MultiQC](https://multiqc.info/)

## Quick Start

1. Clone the repository:
```bash
git clone https://github.com/MasterSwordf1sh/Nextflow-16s-rRNA-analysis
cd Nextflow-16s-rRNA-analysis
```

2. Prepare your input data:
   - Place your paired-end FASTQ files in `data/raw/`
   - Create a metadata file in TSV format and place it at `data/metadata.tsv`

3. Run the pipeline:
```bash
nextflow run main.nf -profile docker
```

## Input

The pipeline requires:
- Paired-end FASTQ files (gzipped recommended)
- A metadata file in TSV format with sample information
- Reference databases (SILVA recommended)

## Output

Results are organized in the following directory structure:

```
results/
├── 1_fastqc/                  # FastQC quality reports
├── 2_trimmed/                 # Trimmed reads
├── 3_qiime2_import/           # QIIME2 imported data
├── 4_dada2/                   # DADA2 denoising results
├── 5_summaries/               # Feature table summaries
├── 6_phylogeny/               # Phylogenetic trees
├── 7_taxonomy/                # Taxonomic classification
├── 8_barplots/                # Taxonomic barplots
├── 9_diversity/               # Diversity metrics and visualizations
├── 10_stats/                  # Statistical analysis results
├── 11_differential_abundance/ # Differential abundance results
└── multiqc_report.html        # MultiQC summary report
```

## Configuration

The pipeline can be configured using the following parameters:

```nextflow
params {
    // Input
    reads = "data/raw/*_R{1,2}.fastq.gz"
    metadata = "data/metadata.tsv"
    
    // Output
    outdir = "results"
    
    // Resources
    cpus = 8
    memory = '16.GB'
    
    // Databases
    silva_db = "databases/silva_138_99_16S.qza"
    
    // Analysis parameters
    // ... (see pipeline script for all options)
}
```

For a complete list of parameters, see the comments in the main pipeline script.

## Customization

You can customize the pipeline execution by providing parameters on the command line:

```bash
nextflow run main.nf --reads 'path/to/reads/*_R{1,2}.fastq.gz' --outdir 'my_results' --cpus 16
```

## Citations

If you use this pipeline, please cite:

- QIIME2: Bolyen E, et al. (2019) Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology, 37: 852–857.
- DADA2: Callahan BJ, et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13: 581-583.
- Nextflow: Di Tommaso P, et al. (2017) Nextflow enables reproducible computational workflows. Nature Biotechnology, 35: 316-319.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

[Contributing](/CONTRIBUTING.md) 

## Contact

If you have any questions or feedback, please open an issue on GitHub or email (mxolisinene@outlook.com) 
