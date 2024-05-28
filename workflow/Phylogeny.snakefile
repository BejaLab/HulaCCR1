
rule pplacer_queries:
    input:
        faa = "metadata/ingroup.faa", # expand("analysis/blastx/{dataset}/{basename}_ingroup.faa", zip, dataset = datasets, basename = datasets.values()),
        exclude = "analysis/phylogeny/selected_ingroup.txt"
    output:
        "analysis/phylogeny/pplacer_queries.fasta"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit grep -vf {input.exclude} -o {output} {input.faa}"

rule add_to_alignment:
    input:
        ref = "analysis/phylogeny/sequences.mafft",
        query = "analysis/phylogeny/pplacer_queries.fasta"
    output:
        "analysis/phylogeny/pplacer_input.mafft"
    conda:
        "envs/mafft.yaml"
    threads:
        20
    shell:
        "mafft --auto --addfragments {input.query} --keeplength --thread {threads} --reorder --auto {input.ref} > {output}"

rule trim_pplace_input:
    input:
        dssp = "analysis/pdb/{acc}.dssp".format(acc = pdb_ref),
        aln = "analysis/phylogeny/pplacer_input.mafft"
    output:
        "analysis/phylogeny/pplacer_input.fasta"
    params:
        min_len = 60
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/trim_tms.py"

rule taxit:
    input:
        tree = "analysis/phylogeny/RAxML_bipartitions.txt",
        info = "analysis/phylogeny/RAxML_info.txt",
        aln = "analysis/phylogeny/sequences.trim"
    output:
        directory("analysis/phylogeny/RAxML.refpkg")
    conda:
        "envs/pplacer.yaml"
    shell:
        "taxit create -l locus_tag -P {output} --tree-file {input.tree} --aln-fasta {input.aln} --tree-stats {input.info}"

rule pplacer:
    input:
        refpkg = "analysis/phylogeny/RAxML.refpkg",
        fasta  = "analysis/phylogeny/pplacer_input.fasta"
    output:
        "analysis/phylogeny/pplacer.jplace"
    conda:
        "envs/pplacer.yaml"
    shell:
        "pplacer -o {output} -c {input.refpkg} {input.fasta}"

rule gappa:
    input:
        "analysis/phylogeny/pplacer.jplace"
    output:
        "output/pplacer.newick"
    conda:
        "envs/gappa.yaml"
    shell:
        "gappa examine graft --jplace-path {input} --fully-resolve --out-dir $(dirname {output})"

rule clade_dist:
    input:
        expand("analysis/blastx/{dataset}/{basename}_filtered.tsv", zip, dataset = datasets, basename = datasets.values())
    output:
        "output/clade_dist.csv"
    params:
        datasets = datasets.keys()
    conda:
        "envs/r.yaml"
    script:
        "scripts/clade_dist.R"

rule phylogeny_sequences:
    input:
        refs = "metadata/BCCRs.faa",
        ingroup = "metadata/selected_ingroup.faa"
    output:
        "analysis/phylogeny/sequences.faa"
    params:
        outgroup = 'StramhfCCR',
        ingroup = 'mgCCR1'
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit grep -nrp {params.outgroup} {input.refs} | seqkit grep -nvrp {params.ingroup} | seqkit seq - {input.ingroup} > {output}"

rule phylogeny_align:
    input:
        "analysis/phylogeny/sequences.faa"
    output:
        "analysis/phylogeny/sequences.mafft"
    conda:
        "envs/mafft.yaml"
    threads:
        20
    shell:
        "mafft --thread {threads} --reorder --auto {input} > {output}"

rule selected_ingroup:
    input:
        "metadata/selected_ingroup.faa"
    output:
        "analysis/phylogeny/selected_ingroup.txt"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -ni {input} | cut -f1 -d@ | tr + '\\n' > {output}"

rule phylogeny_tms:
    input:
        dssp = "analysis/pdb/{acc}.dssp".format(acc = pdb_ref),
        aln = "analysis/phylogeny/sequences.mafft"
    output:
        "analysis/phylogeny/sequences.trim"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/trim_tms.py"

rule phylogeny_raxml:
    input:
        "analysis/phylogeny/sequences.trim"
    output:
        "analysis/phylogeny/RAxML_info.txt",
        "analysis/phylogeny/RAxML_bipartitions.txt"
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
