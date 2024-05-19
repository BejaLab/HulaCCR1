from collections import defaultdict
from pandas import read_csv, DataFrame

input_file = str(snakemake.input)
output_file = str(snakemake.output)

col_names = [ 'record_type', 'cluster', 'seq_length', 'pct_id', 'strand', 'ignore1', 'ignore2', 'cigar', 'gene', 'target' ]
df = read_csv(input_file, delimiter = '\t', names = col_names).drop_duplicates(subset = ['cluster', 'gene'], keep = 'first')
df['gene'] = df['gene'].str.split(' ').str[0]
df['dataset'] = df['gene'].str.split('|').str[0]

df_wide = df.pivot_table(index='cluster', columns='dataset', values='gene', aggfunc=lambda x: ','.join(x) if x.any() else '*').fillna('*').reset_index()

# Count the number of datasets and genes per cluster
num_datasets = df.groupby('cluster')['dataset'].nunique()
num_genes = df.groupby('cluster')['gene'].apply(lambda x: ','.join(x).count(',') + 1 if x.any() else 0)

# Merge the number of datasets and genes back into the wide dataframe
df_wide = df_wide.merge(num_datasets, on='cluster', how='left').rename(columns={'dataset': '# Species'})
df_wide = df_wide.merge(num_genes, on='cluster', how='left').rename(columns={'gene': 'Genes'})

# Sort the dataframe by decreasing num_datasets and num_genes
df_wide = df_wide.drop('cluster', axis = 1).sort_values(by = ['# Species', 'Genes'], ascending = False)
df_wide['Alg.-Conn.'] = -1

first_cols = [ '# Species', 'Genes', 'Alg.-Conn.' ]
df_wide = df_wide[first_cols + [ col for col in df_wide if col not in first_cols ]]

df_wide.to_csv(output_file, sep = "\t", index = False)
