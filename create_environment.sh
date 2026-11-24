#!/bin/bash

# Make sure conda is visible
source ~/miniconda3/etc/profile.d/conda.sh

# Create the conda environment using YAML file
if [ -d "~/miniconda3/envs/microbiome_workbench" ]
then
  echo "Environment for microbiome workbench already exists..."
  echo "Skipping conda env create step."
else
  conda env create -f microbiome_workbench.yml
fi

# Activate environment
conda activate microbiome_workbench

# Download genome reference files

mkdir -p genome_references

## Human genome
if [ -f "genome_references/GCF_000001405.40_GRCh38.p14_genomic.fna" ]
then
  echo "Most up-to-date human genome reference already exists..."
  echo "Skipping human genome reference download."
else
  wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
  gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
  mv GCF_000001405.40_GRCh38.p14_genomic.fna genome_references/
fi

## Mouse genome
if [ -f "genome_references/GCF_000001635.27_GRCm39_genomic.fna" ]
then
  echo "Most up-to-date mouse genome reference already exists..."
  echo "Skipping mouse genome reference download."
else
  wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
  gunzip GCF_000001635.27_GRCm39_genomic.fna.gz
  mv GCF_000001635.27_GRCm39_genomic.fna genome_references/
fi

# Download databases for bioinformatic workflows

mkdir -p  databases

## MetaPhlAn CHOCOPhlAn SGB databases
if [ -d "databases/mpa_vJan26_CHOCOPhlAnSGB_202605" ]
then
  echo "mpa_vJan26_CHOCOPhlAnSGB_202605 database already exists..."
  echo "Skipping mpa_vJan26_CHOCOPhlAnSGB_202605 download."
else
  mkdir -p databases/mpa_vJan26_CHOCOPhlAnSGB_202605
  metaphlan \
    --install \
    --index mpa_vJan26_CHOCOPhlAnSGB_202605 \
    --db_dir databases/mpa_vJan26_CHOCOPhlAnSGB_202605
fi

if [ -d "databases/mpa_vOct22_CHOCOPhlAnSGB_202403" ]
then
  echo "mpa_vOct22_CHOCOPhlAnSGB_202403 database already exists..."
  echo "Skipping mpa_vOct22_CHOCOPhlAnSGB_202403 download."
else
  mkdir -p databases/mpa_vOct22_CHOCOPhlAnSGB_202403
  metaphlan \
    --install \
    --index mpa_vOct22_CHOCOPhlAnSGB_202403 \
    --db_dir databases/mpa_vOct22_CHOCOPhlAnSGB_202403
fi

## HUMAnN databases
if [ -d "databases/humann_databases" ]
then
  echo "HUMAnN database folder already exists..."
  echo "Checking to see if all databases have been downloaded...."
else
  mkdir -p databases/humann_databases
fi

if [ -d "databases/humann_databases/chocophlan" ]
then
  echo "HUMAnN chocophlan database already exists..."
  echo "Skipping HUMAnN chocophlan database download."
else
  humann_databases \
    --download chocophlan full \
    databases/humann_databases/
fi

if [ -d "databases/humann_databases/uniref" ]
then
  echo "HUMAnN uniref database already exists..."
  echo "Skipping HUMAnN uniref database download."
else
  humann_databases \
    --download uniref uniref90_ec_filtered_diamond \
    databases/humann_databases/
fi

if [ -d "databases/humann_databases/utility_mapping" ]
then
  echo "HUMAnN utility mapping database already exists..."
  echo "Skipping HUMAnN utility mapping database download."
else
  humann_databases \
    --download utility_mapping full \
    databases/humann_databases/
fi

# Signal end of setup
echo "*** Microbiome workbench environment setup has completed ***"
