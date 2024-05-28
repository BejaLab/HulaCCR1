from Bio import SeqIO

aln_file = snakemake.input['aln']
dssp_file = snakemake.input['dssp']
output_file = str(snakemake.output)
params = dict(snakemake.params)
min_len = params['min_len'] if 'min_len' in params else 0

helices = [ 'G', 'H', 'I' ]
dssp_seq = ''
helix_start = helix_end = None
with open(dssp_file) as file:
    for line in file:
        if line.lstrip().startswith('#'):
            break
    for pos, line in enumerate(file):
        res = line[13:14]
        ss = line[16:17]
        if ss in helices:
            if helix_start is None:
                helix_start = pos
            helix_end = pos
        dssp_seq += res
assert helix_start is not None, "dssp file malformed"

records = []
ref_seq = ref_pos = None
for record in SeqIO.parse(aln_file, 'fasta'):
    records.append(record)
    if ref_seq is None:
        seq = record.seq.replace('-', '').upper()
        if dssp_seq in seq:
            ref_seq = record.seq
            ref_pos = -seq.index(dssp_seq)

assert ref_seq is not None, f"Reference record not found: '{dssp_seq}'"

aln_start = aln_end = None
for aln_pos, res in enumerate(ref_seq):
    if res != '-':
        if ref_pos == helix_start:
            aln_start = aln_pos
        if ref_pos == helix_end:
            aln_end = aln_pos
        ref_pos += 1

with open(output_file, 'w') as file:
    for record in records:
        record.seq = record.seq[aln_start:aln_end+1]
        ungapped = record.seq.replace('-', '')
        if len(ungapped) >= min_len:
            SeqIO.write(record, file, 'fasta')
