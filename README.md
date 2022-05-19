# bioinfo-variety-scripts
A collection of useful bioinformatics scripts I've written over time.


## correct_sample_file.py

**Description:** Take a sample file and shuffle the first two columns, so that the sample ID is first and the family ID is second. This is useful for [Hail](https://hail.is/).

**Language:** Python

**Usage:** correct_sample_file.py <sample_file_name> <output_file_name> 

## GENs_to_Sorted_8bit_BGENs.sh

**Description:** Take GEN files (eg IMPUTE2 output) and convert them to 8-bit sorted BGEN files organized by chr. Requires [qctool](https://www.well.ox.ac.uk/~gav/qctool_v2/).

**Language:** Bash

**Usage:**  ``` ./GENs_to_Sorted_8bit_BGENs.sh ```

NOTE: The script has a variable where the path to your GEN files should go. It should probably be an imput :/

