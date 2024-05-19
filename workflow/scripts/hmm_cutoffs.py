
hmm_file = str(snakemake.input)
txt_file = str(snakemake.output)
cut = snakemake.params['cut']

with open(hmm_file) as hmm, open(txt_file, 'w') as txt:
    for line in hmm:
        key, *values = line.split()
        if key == "NAME":
            name = values[0]
        elif key == cut:
            cutoff = values[0]
            txt.write(f"{name}\t{cutoff}\n")

