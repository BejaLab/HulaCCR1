
rule blastx_all:
    input:
        expand("analysis/blastx/{dataset}/{basename}_ingroup.fna", zip, dataset = datasets, basename = datasets.values())

rule dload_pdb:
    output:
        "analysis/pdb/{acc}.pdb"
    shell:
        "wget -O {output} https://files.rcsb.org/download/{wildcards.acc}.pdb"

rule dssp:
    input:
        "analysis/pdb/{acc}.pdb"
    output:
        "analysis/pdb/{acc}.dssp"
    conda:
        "envs/dssp.yaml"
    shell:
        "mkdssp {input} > {output}"

rule align_refs:
    input:
        "metadata/BCCRs.faa"
    output:
        "analysis/proteins/BCCRs.mafft"
    conda:
        "envs/mafft.yaml"
    threads:
        20
    shell:
        "mafft --thread {threads} --reorder --auto {input} > {output}"

rule trim_tms:
    input:
        dssp = "analysis/pdb/{acc}.dssp".format(acc = pdb_ref),
        aln = "analysis/proteins/BCCRs.mafft"
    output:
        "analysis/proteins/BCCRs_TM.trim"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/trim_tms.py"

rule trim_tms_strip_gaps:
    input:
        "analysis/proteins/BCCRs_TM.trim"
    output:
        "analysis/proteins/BCCRs_TM.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -g -o {output} {input}"

rule makeblastdb_prot:
    input:
        "{prefix}.faa"
    output:
        "{prefix}.faa.pdb"
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

rule blastx:
    input:
        query = "datasets/{dataset}/{basename}.fna",
        db = "analysis/proteins/BCCRs_TM.faa",
        pdb = "analysis/proteins/BCCRs_TM.faa.pdb"
    output:
        "analysis/blastx/{dataset}/{basename}.txt"
    params:
        evalue = 1e-5
    conda:
        "envs/blast.yaml"
    threads:
        8
    shell:
        "blastx -query {input.query} -db {input.db} -outfmt 6 -out {output} -evalue {params.evalue} -num_threads {threads} -ungapped -comp_based_stats F"

rule filter_blastx:
    input:
        fna = "datasets/{dataset}/{basename}.fna",
        tsv = "analysis/blastx/{dataset}/{basename}.txt"
    output:
        fna = "analysis/blastx/{dataset}/{basename}_filtered.fna",
        tsv = "analysis/blastx/{dataset}/{basename}_filtered.tsv"
    params:
        min_len = 25,
        min_ident = 65
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/filter_blastx.py"

rule filter_blastx_ingroup:
    input:
        "analysis/blastx/{dataset}/{basename}_filtered.fna"
    output:
        "analysis/blastx/{dataset}/{basename}_ingroup.fna"
    params:
        ingroup = "mgCCR1"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit grep -rnp {params.ingroup} {input} > {output}"
