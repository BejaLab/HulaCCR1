from collections import defaultdict
from copy import deepcopy
from Bio import SeqIO
from pandas import read_csv
from copy import deepcopy

faa_file = snakemake.input['faa']
fna_file = snakemake.input['fna']
tblout_file = snakemake.input['tblout']
cutoff_file = snakemake.input['cutoffs']

dataset = snakemake.wildcards['dataset']

out_faa_files = snakemake.output['faa']
out_fna_files = snakemake.output['fna']

markers = snakemake.params['markers']
cutoff_factor = snakemake.params['cutoff_factor']

cutoffs_df = read_csv(cutoff_file, delimiter = '\t', names = [ 'marker', 'cutoff' ])
cutoffs = dict(zip(cutoffs_df.marker, cutoffs_df.cutoff))

def read_tblout(fh):
    for line in fh:
        if not line.startswith('#'):
            target, t_accession, query, q_accession, evalue, score, bias, evalue1, score1, bias1, exp, reg, clu, ov, env, dom, rep, inc, *description = line.split(maxsplit = 18)
            yield target, query, float(evalue), float(score), ' '.join(description)

included_nucl = defaultdict(dict)
included_prot = defaultdict(dict)
with open(tblout_file) as file:
    for prot_target, marker, evalue, score, description in read_tblout(file):
        nucl_target = prot_target.rsplit('_', 1)[0]
        if marker not in included_nucl[nucl_target] and score >= cutoffs[marker] * cutoff_factor:
            included_nucl[nucl_target][marker] = prot_target
            included_prot[prot_target][marker] = nucl_target

faa_fds = {}
fna_fds = {}
for i, marker in enumerate(markers):
    faa_fds[marker] = open(out_faa_files[i], 'w')
    fna_fds[marker] = open(out_fna_files[i], 'w')

for record in SeqIO.parse(faa_file, 'fasta'):
    if record.id in included_prot:
        for marker, nucl_target in included_prot[record.id].items():
            new_record = deepcopy(record)
            new_record.id = f"{dataset}|{nucl_target}"
            SeqIO.write(new_record, faa_fds[marker], 'fasta')
for record in SeqIO.parse(fna_file, 'fasta'):
    if record.id in included_nucl:
        for marker, prot_target in included_nucl[record.id].items():
            new_record = deepcopy(record)
            new_record.id = f"{dataset}|{record.id}"
            SeqIO.write(new_record, fna_fds[marker], 'fasta')

for marker in markers:
    faa_fds[marker].close()
    fna_fds[marker].close()
