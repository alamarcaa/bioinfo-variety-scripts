import pandas as pd
import sys

# This script takes a table of GWAS hits that has been built out of several GWAS hit tables and compiles the results
# so that each variant only shows up once, noting every GWAS hit table that it came from.

args = sys.argv[1:]
if len(args) != 2:
    print("Usage: python generate_unique_hits_table.py <all hits table file> <output_dir>")
    sys.exit(1)

path_to_table = args[0]
path_to_output = args[1]
if not path_to_output.endswith("/"):
    path_to_output += "/"

hits_df = pd.read_csv(path_to_table, sep='\t')

# group rows by hits_df.rsid and keep the row with the lowest value of hits_df.pval
unique_hits_df = hits_df.groupby('rsid').apply(lambda x: x.loc[x.pval.idxmin()])
unique_hits_df = unique_hits_df.sort_values(by=['chromosome', 'position'])


def get_phenotypes(row):
    """ Get a list of phenotypes for a variant (returned as a string), built from all the 'phenotypes' entries that the
    same variant had across the tables """

    hits_phenotype = row.Phenotype
    pheno_list = list(hits_df.loc[hits_df.rsid == row.rsid, 'Phenotype'])
    # make hits_phenotype, the phenotype where the SNP has the lowest p-value, the first in the list
    pheno_list.remove(hits_phenotype)
    pheno_list.insert(0, hits_phenotype)

    pheno_string = ','.join(pheno_list)
    return pheno_string


# Add a 'Phenotype' column to the table, containing the list of phenotypes for each variant
unique_hits_df['Phenotype'] = unique_hits_df.apply(get_phenotypes, axis=1)


def get_genes(row):
    # Get a list of genes for a variant, built from all the 'genes' entries that the same variant had across the tables
    hits_genes_og_list = row.genes.split(',')
    genes_col_content_list = list(hits_df.loc[hits_df.rsid == row.rsid, 'genes'])
    # split each string in gene_list into a list of genes
    gene_list = []
    for string in genes_col_content_list:
        sub_gene_list = string.split(',')
        gene_list.extend(sub_gene_list)
    gene_list = list(set(gene_list))

    for gene in hits_genes_og_list:  # Put the genes that come from the lead variant in front
        gene_list.remove(gene)
        gene_list.insert(0, gene)
    gene_string = ','.join(gene_list)
    return gene_string


# Add a 'genes' column to the table, containing a list of possible target genes for each variant
unique_hits_df['genes'] = unique_hits_df.apply(get_genes, axis=1)

unique_hits_df.to_csv(path_to_output + 'unique_hits_table.tsv', sep='\t', index=False)
