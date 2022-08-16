# bioinfo-variety-scripts
A collection of useful bioinformatics scripts I've written over time.


## correct_sample_file.py

**Description:** Take a sample file and swap the first two columns, so that the sample ID is first and the family ID is
second. This is useful for [Hail](https://hail.is/).

**Language:** Python

**Usage:** ```correct_sample_file.py <sample_file_name> <output_file_name> ```

## GENs_to_Sorted_8bit_BGENs.sh

**Description:** Take GEN files (eg IMPUTE2 output) and convert them to 8-bit sorted BGEN files organized by chromosome.
Requires [qctool](https://www.well.ox.ac.uk/~gav/qctool_v2/).

**Language:** Bash

**Usage:**  ``` ./GENs_to_Sorted_8bit_BGENs.sh ```

NOTE: The script has an internal variable where the path to your GEN files should go. It should probably be an input :/

## plink_pca_from_bgens.sh

**Description:** Take the path to a directory containing .bgen files and do PCA with that data.

**Language:** Bash

**Usage:** # ```./plink_pca_from_bgens.sh <./BGEN_file_folders> <sample_file_name>```

## combine_manhattans.py

**Description:** Combine all the GWAS output files in a directory into a single GWAS output file

**Language:** Python

**Usage:** ```combine_manhattans.py <GWAS data directory> <project_prefix> <steps to run> <p-value threshold>```