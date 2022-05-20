import pandas as pd
import os
import sys
from liftover import get_lifter

# This script creates and fills folders 'bed_files_for_UCSC_browser' and 'all_significant_variant_tables' with .bed
# files and .tsv tables respectively. These are later used to find the closest gene for GWAS hits.
#
# Inputs: 1. path to the folder containing the GWAS hits (e.g. 'GWAS_hits/')
#         2. path to output folder (e.g. 'output/')
# Outputs: 1. .bed files for UCSC browser (e.g. 'output/bed_files_for_UCSC_browser/')
#          2. .tsv tables with all significant variants (e.g. 'output/all_significant_variant_tables/')

args = sys.argv[1:]
if len(args) != 3:
    print("Usage: python generate_bed_files_and_tables.py <path to GWAS hits> <project_name> <path to output folder>")
    sys.exit(1)

files_filepath = args[0]
project_name = args[1]
output_filepath = args[2]


def make_bed_file_hg38(df, label: str, out=os.getcwd()):
    """Function that makes bed file out of GWAS results DataFrame, lifting over results to hg38
    """
    filename = out + '/' + label + "_GWAS_hits_hg38.bed"
    converter = get_lifter('hg19', 'hg38')

    with open(filename, 'w') as out_file:
        for i in df.index:

            chromosome_hg37 = 'chr' + str(df.loc[[i]].chromosome.values[0])
            loc_hg37 = int(df.loc[[i]].position.values[0])
            try:
                liftover_tuple_hg38 = converter.convert_coordinate(chromosome_hg37, loc_hg37)[0]
            except IndexError:
                continue
            chromosome = liftover_tuple_hg38[0]
            loc = str(liftover_tuple_hg38[1])

            DeCODE_rsid = df.loc[[i]].rsid.values[0]
            if str(DeCODE_rsid).startswith('rs'):
                rsid = str(DeCODE_rsid).split(":")[0]
            else:
                rsid = str(DeCODE_rsid)
            out_file.write(chromosome + '\t' + loc + '\t' + '\t' + loc + '\t' + rsid + '\n')


# Create folders if they do not exist
if not os.path.exists(output_filepath):
    os.makedirs(output_filepath)
if not os.path.exists(output_filepath + '/bed_files_for_UCSC_browser'):
    os.mkdir(output_filepath + '/bed_files_for_UCSC_browser')
if not os.path.exists(output_filepath + '/all_significant_variant_tables'):
    os.mkdir(output_filepath + '/all_significant_variant_tables')

for i, file in enumerate(os.listdir(files_filepath)):
    if file.endswith('.txt'):
        print(f"File {i+1} of {len(os.listdir(files_filepath))}")
        f = pd.read_csv(files_filepath + '/' + file, sep="\t", chunksize=100000, index_col=False)
        list_of_downsampled_chunks = []
        for chunk in f:
            list_of_downsampled_chunks.append(chunk[chunk['pval'] < 5.0e-8].copy())
        phenotype = ''.join(file.split('.')[0].split('_')[3:])  # Take phenotype name from file name
        # print(f'#{i}: {phenotype}')
        concat_df = pd.concat(list_of_downsampled_chunks)
        concat_df = concat_df.assign(Phenotype=phenotype)
        concat_df.insert(1, 'rsid_short',
                         [rsid.split(':')[0] if rsid.startswith('rs') else rsid for rsid in concat_df.rsid.to_list()])
        # concat_df = concat_df[concat_df['pval'] <= 5e-8]  # additional pvalue filter
        chrs_with_hits = list(
            concat_df['chromosome'].unique())  # list of chromosomes with p-values that cross threshold
        chrs_with_hits.sort()
        if chrs_with_hits:
            n = concat_df.iloc[0, 9]  # get N from dataframe
            for chromosome in chrs_with_hits:
                # Create bed file for UCSC browser
                make_bed_file_hg38(concat_df[concat_df.chromosome == chromosome].copy(),
                                   label=f'{project_name}_chr{str(chromosome)}_{phenotype}_n{n}_w11_2022',
                                   out=output_filepath + '/bed_files_for_UCSC_browser')
                # Create table of variants that crossed significance threshold ('hit table')
                concat_df[concat_df.chromosome == chromosome].to_csv(output_filepath
                                                                     + f'/all_significant_variant_tables/{project_name}_chr{str(chromosome)}_{phenotype}_n{n}_w11_2022_full_table.tsv',
                                                                     sep='\t', index=False)
