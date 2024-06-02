from Bio import SeqIO
from collections import defaultdict

clstr_file = snakemake.input['clstr']
fasta_file = snakemake.input['fasta']
output_file = str(snakemake.output)

clusters = defaultdict(dict)
with open(clstr_file) as file:
    for line in file:
        if line.startswith('>'):
            prefix, cluster = line.split()
        else:
            num, seq_len, seq_name, *identity = line.split()
            seq_name = seq_name[1:][:-3]
            seq_len = int(seq_len[:-3])
            is_not_xx = 'XX' not in seq_name
            if identity[0] == "*" or is_not_xx:
                clusters[cluster][seq_name] = is_not_xx, seq_len

reps = {}
for cluster, members in clusters.items():
    rep = sorted(members, key = lambda seq_name: members[seq_name], reverse = True)[0]
    reps[rep] = True

with open(output_file, 'w') as file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id in reps:
            SeqIO.write(record, file, 'fasta')
            del reps[record.id]

assert len(reps) == 0, f"Some sequence names were not found: {reps.keys()}"
