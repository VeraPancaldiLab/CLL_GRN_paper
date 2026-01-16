# Time-series RNA-seq of Chronic Lymphocytic Leukaemia (CLL)

## Dataset title
Time-series RNA-seq analysis of Chronic Lymphocytic Leukaemia (CLL) autologous
and monoculture in vitro models

## Authors
Hugo Chenel  
Inserm / Toulouse Cancer Research Center (CRCT)

## Year
2025

## Overview
This dataset contains bulk RNA-seq profiles generated from primary chronic
lymphocytic leukaemia (CLL) samples cultured in vitro under two experimental
conditions:

1. Autologous cultures (AUTO): whole peripheral blood mononuclear cells (PBMCs).
2. Monoculture (B): purified CLL B cells obtained by negative selection.

These cultures are coming from the same patients.
The goal of the study is to characterize transcriptional dynamics and survival-
associated pathways during in vitro CLL progression over a 14-day time course.

## Experimental design
- Organism: Homo sapiens
- Patients: 3 (DAR, SAS, PRG)
- Time points: Day 1, Day 4, Day 8, Day 11, Day 14
- Conditions: Autologous PBMCs and purified CLL B cells
- Replicates: 2 biological replicates per condition, patient and time point
- Total RNA-seq libraries: 60

## Sequencing
- Platform: Illumina NovaSeq 6000
- Flow cell: S1
- Layout: Paired-end
- Read length: 2 × 75 bp
- Library preparation: Illumina Stranded mRNA Prep with IDT UD indexes

## Data availability

### GEO (processed data + metadata)
The dataset is publicly available on NCBI Gene Expression Omnibus (GEO):

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313098

GEO Series accession: **GSE313098**

### SRA (raw sequencing data)
Raw FASTQ files are available via the Sequence Read Archive (SRA):

- BioProject: PRJNA1378255
- Linked from the GEO page via the SRA Run Selector

## Repository structure
- `Data/Raw/`  
  Raw sequencing files (FASTQ.gz). Final archival copies are hosted on SRA.

- `Data/Processed/`  
  Gene-level RNA-seq count matrices (RSEM expected counts).

- `Data/Metadata/`  
  Structured sample metadata and sample informations.

- `Data/Docs/`  
  Human-readable documentation describing the dataset, metadata, and file
  organization.

## Citation
Please cite the dataset using the GEO Series accession:

GSE313098

If a manuscript or preprint is published, the GEO record will be updated with
the corresponding reference.
