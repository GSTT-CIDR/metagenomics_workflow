# CIDR metagenomics workflow
### v3.7.4 - Metagenomics Network Release

## Description
A Snakemake workflow producing a PDF summary of taxonomic classification (centrifuge) and AMR screening. This repo is exclusively for deployment within the metagenomics container environment. 

## Installation requirements
This package has been configured to operate from within a container. Key elements required to work are the bound/mounted data directories.
1. {metagenomics_root}:/mnt
2. {MinKNOW data direcotry}:/data
3. {AWS config dir}:/root/.aws

The metagenomics_root directory contains the databases and configuration files for the pipeline to operate. This includes:
1. A centrifuge index
2. A BLASTn nt database
3. RefSeq database for mapping
4. AMR databases
5. Configuration files for workflow
6. Directories for workflow intermediate outputs and reports

**See directory tree for full description of the library environment**