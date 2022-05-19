import pandas as pd
import os
import sys

# Generates a file "all_hits_table.tsv" containing every SNP from the clustered .bed files.

args = sys.argv[1:]
if len(args) != 4:
    print("Usage: python generate_all_hits_table.py <clustered_bed_dir> <closest_gene_dir> <hit_table_dir> <output_dir>")
    sys.exit(1)

clustered_bed_files_path = args[0]
if not clustered_bed_files_path.endswith("/"):
    clustered_bed_files_path += "/"

bed_files_with_closest_genes_path = args[1]
if not bed_files_with_closest_genes_path.endswith("/"):
    bed_files_with_closest_genes_path += "/"

all_significant_variant_tables_path = args[2]
if not all_significant_variant_tables_path.endswith("/"):
    all_significant_variant_tables_path += "/"

output_path = args[3]


out_df = pd.DataFrame()

for file in os.listdir(clustered_bed_files_path):
    # if the file ends in .bed, and it is not empty
    if file.endswith(".bed.cluster") and os.stat(clustered_bed_files_path + file).st_size != 0:

        # read in the file
        cluster_df = pd.read_csv(clustered_bed_files_path + file, sep="\t", header=None)

        # read the file in 'bed_files_with_closest_gene' folder with the same name
        genes_file = file.replace('.cluster', '.closest')  # change the ending of file from '.cluster' to '.closest'
        closest_gene_df = pd.read_csv(bed_files_with_closest_genes_path + genes_file, sep="\t", header=None)

        # read the hit table in 'all_significant_variant_tables_path' folder
        tsv_filename = '_'.join(file.split('_')[:8]) + '_full_table.tsv'  # corresponding hit table file name
        # TODO: it feels wrong to trust the file name formatting of the hit table will be consistent.
        full_table_df = pd.read_csv(all_significant_variant_tables_path + tsv_filename, sep="\t", index_col=0)

        # For each cluster, add one line to the final hits table
        for i in list(set(cluster_df[4].to_list())):
            cluster_df_i = cluster_df.loc[cluster_df[4] == i]
            full_table_df_i = full_table_df[full_table_df['rsid_short'].isin(cluster_df_i[3])]
            closest_gene_df_i = closest_gene_df[closest_gene_df[3].isin(cluster_df_i[3])]
            gene_list = list(set(closest_gene_df_i[7]))  # make a list of all unique elements of closest_gene_df_i[7]
            gene_list_str = ','.join(gene_list)
            min_pval_df_i = full_table_df_i[full_table_df_i['pval'] == full_table_df_i['pval'].min()].copy()
            min_pval_df_i = min_pval_df_i.iloc[0].to_frame().T  # take only the first row if there are tied p-values
            min_pval_df_i['genes'] = gene_list_str
            out_df = pd.concat([out_df, min_pval_df_i])

    else:
        continue

# sort out_df by column chromosome, then by column position
if 'chromosome' in out_df.columns:
    out_df = out_df.sort_values(by=['chromosome', 'position'])

    output_filepath = output_path + '/all_hits_table.tsv'
    out_df.to_csv(output_filepath, sep="\t", index=False)
else:
    raise ValueError("The output table does not have a 'chromosome' column.")
