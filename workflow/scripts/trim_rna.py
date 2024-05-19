from Bio import SeqIO

tblout_file = snakemake.input['tblout']
fasta_file = snakemake.input['fasta']
output_file = str(snakemake.output)

min_len = snakemake.params['min_len']

regions = {}
with open(tblout_file) as file:
    for line in file:
        if not line.startswith('#'):
            target, target_acc, query, query_acc, mdl, mdl_from, mdl_to, seq_from, seq_to, strand, trunc, pass_num, gc, bias, score, E_value, inc, *description = line.split()
            if target not in regions:
                if strand == '+':
                    regions[target] = int(seq_from) - 1, int(seq_to)
                else:
                    regions[target] = int(seq_to) - 1, int(seq_from)

with open(output_file, 'w') as output:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id in regions:
            start, stop = regions[record.id]
            record = record[start:stop]
        if len(record.seq) >= min_len:
            SeqIO.write(record, output, 'fasta')
