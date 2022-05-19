#!/bin/bash

# Step 1: Get bed files and tables
python generate_bed_files_and_tables.py ../../output/GWAS/ VB_Immune_GWAS .

# Step 2: Assign closest gene to each variant
./closest_gene.sh

# Step 3: Cluster variants together
python generate_all_hits_table.py ./bed_files_clustered/ ./bed_files_with_closest_gene/ ./all_significant_variant_tables/ .
python generate_unique_hits_table.py ./all_hits_table.tsv .
python get_actual_unique_hits_group_by_genes.py ./unique_hits_table.tsv 11 ./unique_hits_group_by_genes.tsv

