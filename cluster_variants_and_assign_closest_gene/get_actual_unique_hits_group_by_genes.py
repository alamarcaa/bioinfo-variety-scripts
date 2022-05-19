import pandas as pd
import functools
import sys

# This script takes the hit table produced by the 'generate_unique_hits_table.py' script and returns ACTUALLY unique 
# hits. It does not work if the genes shared by several hits don't all share at least one common gene.
# It is also overzealous in that it will group together hits that may actually be separate if they share a common
# closest gene.

# Input:
# pval threshold should be an integer, the value of the exponent of the pvalue threshold.
# If you want a 1E-11 threshold, enter 11.

args = sys.argv[1:]
if len(args) != 2:
    print("Usage: python generate_unique_hits_table.py <path_to_hits_table> <pval_threshold> <output_file>")
    sys.exit(1)

path_to_hits_table = args[0]  # For example: 'unique_hits_table_p11.txt'

output_file = args[2]  # For example: 'unique_hits_table_p11_ACTUALLY_unique_genes.txt'

########################################################################################################################
# Extra threshold on hits, specified as an input

pval_threshold = args[1]
pval_threshold = int(pval_threshold)  # make it an int

intermediate_hits_file = f'/home/antton/Projects/Immune_GWAS/scripts/cluster_variants_and_assign_closest_gene' \
                         f'/unique_hits_table_p{pval_threshold}.txt '
# for line in in the file, if the 'pval' column is more than threshold, drop the line
# if the 'pval' column is less than 1E-11, print the line

with open(path_to_hits_table, 'r') as f:
    with open(intermediate_hits_file, 'w') as f2:
        for line in f:
            # skip the first line
            if line.startswith('rsid'):
                f2.write(line)
                continue

            line = line.strip().split('\t')
            if float(line[6]) > float(f'1E-{pval_threshold}'):
                continue
            else:
                f2.write('\t'.join(line) + '\n')


########################################################################################################################
# Group hits by genes


def grouper(sequence):
    # takes a list of lists, returns common elements among them
    group, members = [], set()

    for item in sequence:
        if group and members.isdisjoint(item):
            # new group, yield and start new
            yield group
            group, members = [], set()
        group.append(item)
        members.update(item)

    yield group


# load path_to_hits_table into a df
hits_df = pd.read_csv(intermediate_hits_file, sep='\t')
gene_col_list = list(hits_df.genes.values)
genes_list = [item.split(',') for item in gene_col_list]  # list of lists with all the gene names

grouped_genes_list = []
for gene_group in grouper(genes_list):
    grouped_genes_list.append(gene_group)

print(len(grouped_genes_list))
print(grouped_genes_list)

new_gene_col_list = []
for x, gene_group in enumerate(grouped_genes_list):
    # get common elements in group
    group_size = len(gene_group)
    common_genes = list(functools.reduce(set.intersection, [set(item) for item in gene_group]))
    if not common_genes:  # if there is no gene that is common across all lists, tag it as -undetermined-
        common_genes = [f'-undetermined_{x}-']
    print(group_size, common_genes)
    for i in range(group_size):
        new_gene_col_list.append(','.join(list(common_genes)))

print(new_gene_col_list)
# group hits_df by grouped_genes_list, keep only the row with the lowest p-value
# and keep the gene name
hits_df['genes'] = new_gene_col_list
hits_df['pval'] = hits_df['pval'].astype(float)  # convert pval to float
print(hits_df.head(30))
# Group by 'genes', keep only the row with the lowest value in the 'pval' column
unique_hits_df = hits_df.loc[hits_df.groupby('genes').pval.idxmin()]
# sort by index
unique_hits_df = unique_hits_df.sort_index()
print(unique_hits_df)

unique_hits_df.to_csv(output_file, sep='\t', index=False)
