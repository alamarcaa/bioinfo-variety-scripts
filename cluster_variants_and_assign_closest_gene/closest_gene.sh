#!/bin/bash
# This script goes over a bunch of .bed files and for each SNP gets the closest gene.
# It requires a folder named 'bed_files_for_UCSC_browser' containing all the .bed files, as well as a reference gene list.
# It outputs two folders: 'bed_files_with_closes_gene' and 'bed_files_clustered', each containing the same amount of
# files as the bed_files_for_UCSC_browser folder.

# Usage: ./closest_gene.sh (This script is meant to be run in the '<project>/data/interim' folder.)

# Inspired by: https://bioinformatics.stackexchange.com/questions/3321/finding-the-nearest-gene-to-a-specified-gene-region

# aliases
bedtools='/home/antton/Programs/bedtools2/bin/bedtools'

# This script requires a reference gene list with the columns: chromosome, start, end, gene_name. This file is expected
# to be named 'NCBI_refseq_all_genes_sorted.txt' and be located in the currect directory.
# The following section takes the NCBI_refseq_all_genes list (you can download it from here:
# https://genome.ucsc.edu/cgi-bin/hgTables)9 and cleans it up so that we take the chrom, txStart, txEnd, name2 columns.
# It then replaces the spaces for tabs, takes away the header, and sorts the file with bedtools.
# Finally, it removes all the intermediary files. If you have already produced the cleaned table, this section is skipped.
if [ ! -f NCBI_refseq_all_genes_sorted.txt ]; then
    awk '{print $3, $5, $6, $13}' ../../raw/NCBI_refseq_all_genes > NCBI_refseq_all_genes_cleaned.txt
    # replace spaces with tabs for file NCBI_refseq_all_genes_cleaned.txt
    sed 's/ /\t/g' NCBI_refseq_all_genes_cleaned.txt > NCBI_refseq_all_genes_cleaned_tab.txt
    # remove header
    tail -n +2 NCBI_refseq_all_genes_cleaned_tab.txt > NCBI_refseq_all_genes_cleaned_tab_noheader.txt
    # sort
    $bedtools sort -i NCBI_refseq_all_genes_cleaned_tab_noheader.txt > NCBI_refseq_all_genes_sorted.txt
    # remove all the intermediary files
    rm NCBI_refseq_all_genes_cleaned.txt NCBI_refseq_all_genes_cleaned_tab.txt NCBI_refseq_all_genes_cleaned_tab_noheader.txt
fi

# If folders do not exist, create them
if [ ! -d "bed_files_interim" ]; then
    mkdir bed_files_interim
fi
if [ ! -d "bed_files_with_closest_gene" ]; then
    mkdir bed_files_with_closest_gene
fi
if [ ! -d "bed_files_clustered" ]; then
    mkdir bed_files_clustered
fi

# For each file in folder bed_files_for_UCSC_browser, replace tabs, sort and intersect with
# NCBI_refseq_all_genes_sorted.txt. Get the closest gene to each SNP, and cluster SNPs within 100000bp.
for file in bed_files_for_UCSC_browser/*.bed
do
  filename=$(basename "$file")
  sed 's/\t\t/\t/g' "$file" > bed_files_interim/$filename.tab
  tail -n +2 bed_files_interim/$filename.tab > bed_files_interim/$filename.tab_noheader
  $bedtools sort -i bed_files_interim/$filename.tab_noheader > bed_files_interim/$filename.sorted
  $bedtools closest -a bed_files_interim/$filename.sorted -b NCBI_refseq_all_genes_sorted.txt > bed_files_with_closest_gene/$filename.closest
  $bedtools cluster -d 100000 -i bed_files_interim/$filename.sorted > bed_files_clustered/$filename.cluster

done

rm -r bed_files_interim