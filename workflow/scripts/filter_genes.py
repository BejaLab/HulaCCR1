from collections import defaultdict
from pandas import read_csv, DataFrame

tsv_file = snakemake.input['tsv']
clades_file = snakemake.input['clades']

output_file = str(snakemake.output)

assemblage = snakemake.params['assemblage']
reference = snakemake.params['reference']
min_datasets = snakemake.params['min_datasets']
all_datasets = snakemake.params['datasets']

datasets = defaultdict(list)
with open(clades_file) as file:
    for line in file:
        dataset, clade, num = line.split(',')
        clade = clade.split('=')[0].split('|')[2]
        if clade in assemblage and dataset in all_datasets:
            datasets[dataset].append(clade)
datasets = { dataset: clades for dataset, clades in datasets.items() if len(clades) == len(assemblage) }

df = read_csv(tsv_file, delimiter = '\t', header = 0)

columns = {}
ref_col = None
for col in df.columns:
    for dataset in datasets:
        if col.startswith(dataset):
            columns[col] = dataset
    if col.startswith(reference):
        ref_col = col

genes = []
if ref_col:
    for orthogroup, row in df.iterrows():
        if row[ref_col] != '*':
            non_empty = { col: dataset for col, dataset in columns.items() if row[col] != '*' }
            if len(non_empty) / len(datasets) >= min_datasets:
                for col, dataset in non_empty.items():
                    for gene in row[col].split(','):
                        genes.append({ "gene": gene, "dataset": dataset, "orthogroup": orthogroup, "num_datasets": len(non_empty), "fract_datasets": len(non_empty) / len(datasets) })
output = DataFrame(genes)
output.to_csv(output_file, sep = "\t", index = False)
