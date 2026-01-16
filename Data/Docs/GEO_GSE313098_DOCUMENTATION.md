# GEO Submission Documentation – GSE313098

## Accession
GSE313098

## GEO link
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313098

## Status
- Submission date: 09 Dec 2025
- Last update: 10 Dec 2025
- Scheduled public release: 10 Dec 2025
- Status: Approved

## Study summary
This GEO Series contains time-series bulk RNA-seq data from chronic lymphocytic
leukaemia (CLL) primary samples cultured in vitro under two conditions:
autologous PBMC cultures and purified CLL B-cell monocultures.

Samples were collected from three CLL patients (DAR, SAS, PRG) at five time
points (Day 1, Day 4, Day 8, Day 11, Day 14), with two biological replicates per
time point and condition, resulting in 60 RNA-seq samples.

## Processed data in GEO
The following processed files are provided as supplementary data:

- `GSE313098_Gene_count_matrix_ENS_HGCN_auto.txt.gz`  
  Raw RSEM expected gene counts for autologous PBMC samples

- `GSE313098_Gene_count_matrix_ENS_HGCN_B.txt.gz`  
  Raw RSEM expected gene counts for monoculture CLL B-cell samples

File format:
- Tab-delimited text
- Rows: Ensembl gene IDs with HGNC symbols
- Columns: sample identifiers
- Values: raw, unnormalized RSEM expected counts

These files are available via both FTP and HTTP.

## Raw data
All raw sequencing data are hosted by the Sequence Read Archive (SRA):

- BioProject: `PRJNA1378255`
- FASTQ files available via the SRA Run Selector linked on the GEO page

## Sample metadata
All 60 samples were assigned GSM accessions.

Metadata includes:
- Patient identifier
- Time point
- Culture condition
- Biological replicate
- Sequencing platform
- Library preparation details
- Sample-specific attributes (cell type, tissue, molecule type)

## Contact information
Hugo Chenel  
Evotec  
195 Route d'Espagne, 31100 Toulouse, France  
Email: hugo.chenel@evotec.com

Vera Pancaldi  
Toulouse Cancer Research Center (CRCT), Inserm  
2 Avenue Hubert Curien, 31100 Toulouse, France  
Email: vera.pancaldi@inserm.fr

## Notes
This dataset must be manually released earlier if referenced in a manuscript,
preprint, thesis, or other public document prior to the scheduled GEO release
date.
