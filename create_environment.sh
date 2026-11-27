#!/bin/bash

# Make sure conda is visible
source ~/miniconda3/etc/profile.d/conda.sh

# Create the conda environment using YAML file
conda env create -f microbiome_workbench.yml

# Activate environment
conda activate microbiome_workbench

# Download genome reference files

mkdir -p ~/miniconda3/envs/microbiome_workbench/opt/genome_references

## Human genome
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna \
    ~/miniconda3/envs/microbiome_workbench/opt/genome_references/

## Mouse genome
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
gunzip GCF_000001635.27_GRCm39_genomic.fna.gz
mv GCF_000001635.27_GRCm39_genomic.fna \
    ~/miniconda3/envs/microbiome_workbench/opt/genome_references/

# Download databases for bioinformatic workflows

## MetaPhlAn CHOCOPhlAn SGB databases
metaphlan \
    --install \
    --index mpa_vJan26_CHOCOPhlAnSGB_202605

metaphlan \
    --install \
    --index mpa_vOct22_CHOCOPhlAnSGB_202403

## HUMAnN databases
mkdir -p ~/miniconda3/envs/microbiome_workbench/opt/humann_databases

humann_databases \
    --download chocophlan full \
    ~/miniconda3/envs/microbiome_workbench/opt/humann_databases/

humann_databases \
    --download uniref uniref90_ec_filtered_diamond \
    ~/miniconda3/envs/microbiome_workbench/opt/humann_databases/

humann_databases \
    --download utility_mapping full \
    ~/miniconda3/envs/microbiome_workbench/opt/humann_databases/


echo " "
echo "*** Microbiome workbench environment setup has completed ***"
echo " "