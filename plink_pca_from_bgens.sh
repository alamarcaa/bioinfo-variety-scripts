#!/bin/bash

##Shell script that takes the path to a directory containing .bgen files and does a PCA for those
#Requires: python, cat-bgen, bgenix, plink, as well as the "correct_sample_file.py" script
#INPUT:
# $1 : input path, path to BGENs folder  $2 : sample file
#EXAMPLE:
# ./plink_pca_from_bgens.sh './BGEN_files' 'Sweden_CordBlood_wAliases_211029.sample'

#TODO: output to log instead of screen?
#TODO: IDK if doing the sample file correction here is reasonable:/

#Aliases
catbgen='/home/antton/Programs/bgen/bin/cat-bgen'
bgenix='/home/antton/Programs/bgen/bin/bgenix'
plink='/home/antton/Programs/plink2'

#1 Create 'bgen_list.txt'
  #Check if directory exists
if [ ! -d "$1" ]; then
  echo "No directory '$1' could be found."
  exit
fi

[ -f bgen_list.txt ] && { rm bgen_list.txt; } #If list exists delete it to avoid problems 
for i in {1..22}; do
  [ -f $1/chr$i.bgen ] && { echo $1/chr$i.bgen >> bgen_list.txt; } #Check that each BGEN exists and add its name to file if it does 
done

#2 Concatenate all bgens
mkdir merged_BGEN
$catbgen -g $(cat bgen_list.txt) -og merged_BGEN/merged.bgen #NOTE: my local alias is "cat-bgen" with a '-' !!

#3 Index the merged BGEN file, then produced a sorted file using 'bgenix'
$bgenix -g merged_BGEN/merged.bgen -index
$bgenix -g merged_BGEN/merged.bgen > merged_BGEN/merged_sorted.bgen

#4 Correct sample file
CORRECTED_SAMPLE="CORRECTED_$2"
python correct_sample_file.py "$2" "$CORRECTED_SAMPLE"

#5 Do PCA

mkdir PCA_results
$plink --make-rel --bgen ./merged_BGEN/merged_sorted.bgen --sample "$CORRECTED_SAMPLE" --missing-code -9,NA --maf 0.05 --hwe 1e-6 --indep-pairwise 1000 0.3 --pca 20 --out ./PCA_results/PCA_from_BGENs


