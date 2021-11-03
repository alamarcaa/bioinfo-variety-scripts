#!/usr/bin/bash
#
# Script that takes GEN files outputed by IMPUTE2 and converst them into 8bit sorted BGEN files suitable for HAIL
#
# Author: Antton Lamarca, October 2021 

echo "Commencing... keep in mind the conversion will take some time!"

# Step 0: check we have the right path to the data
gen_filepath="deCODE_IMPUTE2_GEN_files" #Name of the folder where the GEN files are
[ -d $gen_filepath/ ] && echo "Directory $gen_filepath exists."

# Step 1: separate GEN files into folders by chromosome.
# This step also creates txt file with the names of 22 chromosomes at the directory where it's being run.
# Note that the generated subfolders will contain GEN files that are not sorted.

echo "Commencing Step 1: separating GEN files into folders by chromosome."

mkdir GEN_files_by_chromosomes # Directory that will contain all GEN files in sub-folders for each chromosome
cd GEN_files_by_chromosomes
for num in {1..22}; do
  # echo chr$num >> chr_list.txt
  mkdir chr$num
  mv ../$gen_filepath/*chr$num.*.impute2.cl.gz chr$num/
done
cd ..

echo "Step 1 complete! GEN files separated."

# Step 2: per chromosome, fuse all .gz files in a single chr folder. It outputs the file inside of the subfolder itslef.

echo "Commencing Step 2: concatenating GEN files into singe GEN files per chromosome."

mkdir fused_GEN_files # Directory that will contain concatenated GEN files, one per chromosome
cd fused_GEN_files

counter=1
for folder in ../GEN_files_by_chromosomes/*/; do 
    echo "Sorting folder $folder...";
    var=$(ls $folder | sort -V);
    cd $folder;
    filename=$(echo $folder | cut -d"/" -f3); #variable name is path dependant! Be careful with 'cut' command
    echo "Fusing all files into a single $filename.gz file...";
    cat $var > ../../fused_GEN_files/$filename.gz;
    echo $counter/22 complete;
    (( counter++ ))
    cd ../../fused_GEN_files
done

echo "Step 2 complete! GEN files have been concatenated."

# Step 3: make each chromosome GEN file into a sorted 8bit BGEN file

# In order for qctool to be callable, make sure to give it an alias by using:
# $alias qctool='/home/antton/Programs/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool'

echo "Commencing Step 3: converting GEN files into sorted BGEN files using qctool."

cd ..
mkdir BGEN_files # Directory that will contain the 8bit BGEN files.

for file in fused_GEN_files/*; do
    filename=$(echo $file | cut -d"/" -f2);
    chrname="${filename%.*}";
    chrnumber=$(echo $chrname | cut -c 4-); #Relative path to qctool should work fine here I think
    /home/antton/Programs/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool -bgen-bits 8 -assume-chromosome $chrnumber -g fused_GEN_files/$chrname.gz -filetype gen -og BGEN_files/$chrname.bgen -ofiletype bgen -sort
done #Notice the '-sort' flag! QCTOOL sorts based on on genomic position, alleles, and ID fields.


echo "Step 3 complete! BGENs generated."
