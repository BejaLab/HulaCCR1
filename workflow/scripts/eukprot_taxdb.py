from Bio import SeqIO
from pandas import read_csv
import sys

metadata_file = snakemake.input['metadata'] # "EukProt_included_data_sets.v03.2021_11_22.txt"
unieuk_file = snakemake.input['unieuk'] # "unieuk-0.0.1-pre_release.tsv"
fasta_file = snakemake.input['fasta'] # "proteins.fasta"

names_file = snakemake.output['names']
nodes_file = snakemake.output['nodes']
mapping_file = snakemake.output['mapping']

metadata = read_csv(metadata_file, delimiter = '\t').set_index('EukProt_ID').T.to_dict()

def write_node(fh, tax_id, parent_id, rank):
    fields = [ str(tax_id), str(parent_id), rank ]
    fields += '' * 10
    fh.write('\t|\t'.join(fields) + '\t|\n')

def write_name(fh, tax_id, taxon):
    fields = [ str(tax_id), taxon, f"{taxon} [{tax_id}]", 'scientific name' ]
    fh.write('\t|\t'.join(fields) + '\t|\n')

def get_eukprot(eukprot_record):
    rank = eukprot_record['rank'] if eukprot_record['rank'] else 'clade'
    tax_id = unieuk[lineage]['UniEuk id'].split('#')[1]
    return rank, int(tax_id)

tax_ids = {}
last_tax_id = 1000000
eukprot_tax = {}

unieuk = read_csv(unieuk_file, na_filter = False, delimiter = '\t').set_index('taxon').T.to_dict()

with open(names_file, 'w') as names, open(nodes_file, 'w') as nodes:
    write_name(names, 1, 'root')
    write_node(nodes, 1, 1, 'no rank')
    for eukprot_id, data in metadata.items():
        taxa = data['Taxonomy_UniEuk'].split(';')
        parent_id = 1
        lineage = ''
        for taxon in taxa:
            if taxon.startswith('strain'):
                next
            lineage += taxon + ';'
            if lineage in tax_ids:
                tax_id = tax_ids[lineage]
            else:
                rank = ''
                if lineage in unieuk:
                    rank, tax_id = get_eukprot(unieuk[lineage])
                else:
                    last_tax_id += 1
                    tax_id = last_tax_id
                    rank = 'clade'
                tax_ids[lineage] = tax_id
                write_name(names, tax_id, taxon)
                write_node(nodes, tax_id, parent_id, rank)
            parent_id = tax_id
        eukprot_tax[eukprot_id] = parent_id

with open(mapping_file, 'w') as mapping:
    with open(fasta_file) as fasta:
        for line in fasta:
            if line[0] == '>':
                prot_id = line.split()[0][1:]
                eukprot_id = prot_id[:7]
                mapping.write(f"{prot_id}\t{eukprot_tax[eukprot_id]}\n")
