
rule all_flanks:
    input:
        "output/transposons_pplacer.newick",
        "output/transposons_subclades.csv"

rule select_flanks:
    input:
        "metadata/flanks.fna"
    output:
        "analysis/flanks/flanks.fna"
    params:
        min_len = 150
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -gm {params.min_len} -o {output} {input}"

rule cluster_flanks:
    input:
        "analysis/flanks/flanks.fna"
    output:
        "analysis/flanks/flanks.cdhit"
    params:
        c = 0.9
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output} -c {params.c} -d 0"

rule cluster_makeblastdb:
    input:
        "analysis/flanks/flanks.cdhit"
    output:
        "analysis/flanks/flanks.cdhit.ndb"
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule flanks_blastn:
    input:
        db = "analysis/flanks/flanks.cdhit",
        ndb = "analysis/flanks/flanks.cdhit.ndb",
        query = "datasets/{dataset}/{basename}.fna"
    output:
        "analysis/flanks/blastn/{dataset}_{basename}.txt"
    conda:
        "envs/blast.yaml"
    shell:
        "blastn -query {input.query} -db {input.db} -out {output} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'"

rule flanks_blastn_filter:
    input:
        "analysis/flanks/blastn/{dataset}_{basename}.txt"
    output:
        "analysis/flanks/blastn/{dataset}_{basename}_filtered.tsv"
    params:
        evalue = 1e-12
    conda:
        "envs/kits.yaml"
    shell:
        "awk -ve={params.evalue} '$11<=e' {input} > {output}"

rule flanks_blastn_extract:
    input:
        fna = "datasets/{dataset}/{basename}.fna",
        tsv = "analysis/flanks/blastn/{dataset}_{basename}_filtered.tsv"
    output:
        "analysis/flanks/blastn/{dataset}_{basename}.fna"
    conda:
        "envs/kits.yaml"
    shell:
        "cut -f1 {input.tsv} | seqkit grep -f- {input.fna} -o {output}"

rule transposonpsi:
    input:
        "analysis/flanks/blastn/{dataset}_{basename}.fna"
    output:
        "analysis/flanks/tpsi/{dataset}_{basename}.gff3"
    shadow:
        "shallow"
    conda:
        "envs/transposonpsi.yaml"
    shell:
        "transposonPSI.pl {input} nuc && mv $(basename {input}).TPSI.allHits.chains.bestPerLocus.gff3 {output}"

rule transposonpsi_filter:
    input:
        "analysis/flanks/tpsi/{dataset}_{basename}.gff3"
    output:
        "analysis/flanks/tpsi/{dataset}_{basename}_filter.gff"
    params:
        evalue = 1e-10
    conda:
        "envs/kits.yaml"
    shell:
        "awk 'match($0,/E=(\\S+)/,E) && E[1]<{params.evalue}' {input} > {output}"

rule transposonpsi_fna:
    input:
        gff = "analysis/flanks/tpsi/{dataset}_{basename}_filter.gff",
        fna = "analysis/datasets/{dataset}_{basename}.fna",
        fai = "analysis/datasets/{dataset}_{basename}.fna.fai"
    output:
        "analysis/flanks/tpsi/{dataset}_{basename}.fna"
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools getfasta -s -fi {input.fna} -bed {input.gff} > {output}"

rule transposonpsi_faa:
    input:
        "analysis/flanks/tpsi/{dataset}_{basename}.fna"
    output:
        "analysis/flanks/tpsi/{dataset}_{basename}.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit translate {input} | seqkit replace -p ^ -r '{wildcards.dataset}|' -o {output}"

rule blastp_ref_databases:
    input:
        query = "analysis/flanks/tpsi/{dataset}_{basename}.faa",
        db_dir = "databases/{database}"
    output:
        "analysis/flanks/databases/{dataset}_{basename}_{database}.txt"
    params:
        evalue = 1e-10,
        db = lambda w: "databases/{database}/{prefix}".format(database = w.database, prefix = databases[w.database])
    conda:
        "envs/blast.yaml"
    shell:
        "blastp -query {input.query} -db {params.db} -out {output} -outfmt 6 -evalue {params.evalue}"

rule blastp_ref_databases_filter:
    input:
        "analysis/flanks/databases/{dataset}_{basename}_{database}.txt"
    output:
        "analysis/flanks/databases/{dataset}_{basename}_{database}_filtered.tsv"
    params:
        score = 250
    conda:
        "envs/kits.yaml"
    shell:
        "awk '$12>={params.score}' {input} > {output}"

rule blastp_ref_databases_extract:
    input:
        tsv = expand("analysis/flanks/databases/{dataset}_{basename}_{{database}}_filtered.tsv", zip, dataset = datasets, basename = datasets.values()),
        db_dir = "databases/{database}"
    output:
        "analysis/flanks/refs/{database}.faa"
    params:
        db = lambda w: "databases/{database}/{prefix}".format(database = w.database, prefix = databases[w.database])
    conda:
        "envs/blast.yaml"
    shell:
        "cut -f2 {input.tsv} | sort -u | blastdbcmd -db {params.db} -entry_batch - | sed 's/>/>{wildcards.database}|/' > {output}"

rule blastp_ref_databases_cat:
    input:
        expand("analysis/flanks/refs/{database}.faa", database = databases)
    output:
        "analysis/flanks/refs.faa"
    shell:
        "cat {input} > {output}"

rule blastp_ref_databases_cdhit:
    input:
        "analysis/flanks/refs.faa"
    output:
        "analysis/flanks/refs.cdhit"
    conda:
        "envs/cd-hit.yaml"
    params:
        c = 0.62,
        n = 2
    shell:
        "cd-hit -i {input} -o {output} -c {params.c} -n {params.n}"

rule datasets_pol_cat:
    input:
        expand("analysis/flanks/tpsi/{dataset}_{basename}.faa", zip, dataset = datasets, basename = datasets.values())
    output:
        "analysis/flanks/datasets.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "cat {input} > {output}"

rule datasets_pol_cdhit:
    input:
        "analysis/flanks/datasets.faa"
    output:
        fasta = "analysis/flanks/datasets.cdhit",
        clstr = "analysis/flanks/datasets.cdhit.clstr"
    params:
        c = 0.95
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output.fasta} -c {params.c} -d 0"

rule transposons_cat:
    input:
        datasets = "analysis/flanks/datasets.cdhit",
        refs = "analysis/flanks/refs.cdhit"
    output:
        "analysis/flanks/transposons.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit replace -sp [^ACDEFGHIKLMNPQRSTVWY] -r X -o {output} {input}"

rule transposons_mafft:
    input:
        "analysis/flanks/transposons.faa"
    output:
        "analysis/flanks/transposons.mafft"
    conda:
        "envs/mafft.yaml"
    threads:
        10
    shell:
        "mafft --thread {threads} --auto {input} > {output}"

rule transposons_trimal:
    input:
        "analysis/flanks/transposons.mafft"
    output:
        "analysis/flanks/transposons.trimal"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -in {input} -out {output} -automated1 -keepheader"

rule transposons_long:
    input:
        "analysis/flanks/transposons.trimal"
    output:
        "analysis/flanks/transposons.long"
    params:
        m = 200
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -gm {params.m} {input} | seqkit seq -ni | seqkit grep -f- {input} | seqkit replace -p '[:()]' -r _ -o {output}"

rule transposons_raxml:
    input:
        "analysis/flanks/transposons.long"
    output:
        "analysis/flanks/RAxML_info.txt",
        "analysis/flanks/RAxML_bipartitions.txt",
        "analysis/flanks/RAxML_bestTree.txt",
        "analysis/flanks/RAxML_bipartitionsBranchLabels.txt",
        "analysis/flanks/RAxML_bootstrap.txt"
    params:
        model = "PROTCATLG",
        seed = 123,
        bootstrap = 1000
    conda:
        "envs/raxml.yaml"
    threads:
        20
    shell:
        "raxmlHPC-PTHREADS-SSE3 -f a -p {params.seed} -x {params.seed} -# {params.bootstrap} -m {params.model} -T {threads} -s {input} -n txt -w $(dirname $(realpath {output}))"

rule transposons_taxit:
    input:
        tree = "analysis/flanks/RAxML_bipartitions.txt",
        info = "analysis/flanks/RAxML_info.txt",
        aln = "analysis/flanks/transposons.long"
    output:
        directory("analysis/flanks/RAxML.refpkg")
    conda:
        "envs/pplacer.yaml"
    shell:
        "taxit create -l locus_tag -P {output} --tree-file {input.tree} --aln-fasta {input.aln} --tree-stats {input.info}"

rule transposons_pplacer_input:
    input:
        "analysis/flanks/transposons.trimal"
    output:
        "analysis/flanks/transposons.trimal.fasta"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit replace -p '[:()]' -r _ -o {output} {input}"

rule transposons_pplacer:
    input:
        refpkg = "analysis/flanks/RAxML.refpkg",
        fasta  = "analysis/flanks/transposons.trimal.fasta"
    output:
        "analysis/flanks/transposons_pplacer.jplace"
    conda:
        "envs/pplacer.yaml"
    shell:
        "pplacer -o {output} -c {input.refpkg} {input.fasta}"

rule transposons_gappa:
    input:
        "analysis/flanks/transposons_pplacer.jplace"
    output:
        "output/transposons_pplacer.newick"
    conda:
        "envs/gappa.yaml"
    shell:
        "gappa examine graft --jplace-path {input} --fully-resolve --out-dir $(dirname {output})"

rule transposons_subclades:
    input:
        tsv = expand("analysis/flanks/blastn/{dataset}_{basename}_filtered.tsv", zip, dataset = datasets, basename = datasets.values()),
        clstr = "analysis/flanks/datasets.cdhit.clstr"
    output:
        "output/transposons_subclades.csv"
    conda:
        "envs/r.yaml"
    script:
        "scripts/transposons_subclades.R"
