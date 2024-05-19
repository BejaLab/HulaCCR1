from Bio import SeqIO

tsv_input = snakemake.input['tsv']
fna_input = snakemake.input['fna']

tsv_output = snakemake.output['tsv']
fna_output = snakemake.output['fna']

min_len = snakemake.params['min_len']
min_ident = snakemake.params['min_ident']

seen = {}
sseqids = {}
revcomp = {}
with open(tsv_input) as fh_in, open(tsv_output, 'w') as fh_out:
    for line in fh_in:
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
        if qseqid not in seen:
            if int(length) >= min_len and float(pident) >= min_ident:
                sseqids[qseqid] = sseqid
                revcomp[qseqid] = int(qstart) > int(qend)
                fh_out.write(line)
            seen[qseqid] = True

with open(fna_output, 'w') as fh_out:
    for record in SeqIO.parse(fna_input, 'fasta'):
        if record.id in sseqids:
            record.description = sseqids[record.id]
            if revcomp[record.id]:
                record.seq = record.seq.reverse_complement()
            SeqIO.write(record, fh_out, 'fasta')
