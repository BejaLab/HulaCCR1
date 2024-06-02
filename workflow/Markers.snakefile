rule metaxa:
    input:
        "datasets/{dataset}/{basename}.fna"
    output:
        "analysis/markers/rna/{gene}/{dataset}_{basename}.extraction.fasta"
    params:
        prefix = "analysis/markers/rna/{gene}/{dataset}_{basename}"
    conda:
        "envs/metaxa.yaml"
    shell:
        "metaxa2_x -i {input} -o {params.prefix} -g {wildcards.gene} -t e --allow_single_domain 1 --summary 0 --graphical 0 --silent 1"

rule metaxa_rename:
    input:
        "analysis/markers/rna/ssu/{dataset}_{basename}.extraction.fasta"
    output:
        "analysis/markers/rna/ssu/{dataset}_{basename}.fna"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit replace -p '[|]' -r ' ' {input} | seqkit replace -p ^ -r '{wildcards.dataset}|' -o {output}"

rule dload_pr2:
    output:
        "analysis/pr2/{filename}.fasta"
    params:
        url_prefix = "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0"
    shell:
        "wget -O- {params.url_prefix}_$(basename {output}).gz | gzip -cd > {output}"

rule sintax:
    input:
        db = "analysis/pr2/SSU_UTAX.fasta",
        query = "analysis/markers/rna/ssu/{dataset}_{basename}.fna"
    output:
        "analysis/markers/rna/ssu/{dataset}_{basename}.sintax.tsv"
    threads:
        4
    shell:
        "usearch -sintax {input.query} -db {input.db} -tabbedout {output} -strand both -threads {threads}"

ingroup = "Kathablepharidacea"
outgroups = [
    "AB508339.1", # Palpitomonas bilix
    "X57162.1", # Guillardia theta
    "LC647572.1", # Hemiselmis andersenii
    "F952575.1", # Rhodomonas lens
    "LC151286.1", # Hemiarma marina
    "AF508277.1" # Goniomonas pacifica CCMP1869
]

rule ssu_pr2:
    input:
        fna = "analysis/pr2/SSU_UTAX.fasta",
        chim = "metadata/kathablepharid_ssu_pr2_chimeric.txt"
    output:
        "analysis/markers/ssu_phylogeny/pr2.fna"
    params:
        incl = '|'.join([ ingroup ] + outgroups),
        excl = '|'.join([ g + ':' for g in outgroups ]),
        min_len = 1400
    conda:
        "envs/kits.yaml"
    shell:
        """
            seqkit grep -rp '{params.incl}' {input.fna} | seqkit grep -vrp '{params.excl}' | \
                seqkit grep -vrf {input.chim} | \
                seqkit seq -gm{params.min_len} | \
                seqkit replace -p 'tax=.+,s:' | \
                seqkit replace -p ';' -r '_' -o {output}
        """

rule ssu_pr2_cluster:
    input:
        "analysis/markers/ssu_phylogeny/pr2.fna"
    output:
        cdhit = "analysis/markers/ssu_phylogeny/pr2.cdhit",
        clstr = "analysis/markers/ssu_phylogeny/pr2.cdhit.clstr"
    params:
        c = 0.99
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output.cdhit} -c {params.c} -d 0"

rule ssu_pr2_cluster_reps:
    input:
        clstr = "analysis/markers/ssu_phylogeny/pr2.cdhit.clstr",
        fasta = "analysis/markers/ssu_phylogeny/pr2.fna"
    output:
        "analysis/markers/ssu_phylogeny/pr2.cdhit.reps"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/select_ssu_representatives.py"

rule ssu_raxml:
    input:
        "analysis/markers/ssu_phylogeny/ssu.trim"
    output:
        "analysis/markers/ssu_phylogeny/RAxML_info.txt",
        "analysis/markers/ssu_phylogeny/RAxML_bipartitions.txt",
        "analysis/markers/ssu_phylogeny/RAxML_bestTree.txt",
        "analysis/markers/ssu_phylogeny/RAxML_bipartitionsBranchLabels.txt",
        "analysis/markers/ssu_phylogeny/RAxML_bootstrap.txt"
    params:
        model = "GTRCAT",
        seed = 123,
        bootstrap = 1000
    conda:
        "envs/raxml.yaml"
    threads:
        20
    shell:
        "raxmlHPC-PTHREADS-SSE3 -f a -p {params.seed} -x {params.seed} -# {params.bootstrap} -m {params.model} -T {threads} -s {input} -n txt -w $(dirname $(realpath {output}))"

rule ssu_align:
    input:
        "metadata/kathablepharid_ssu_Ga0334989.fna",
        "metadata/kathablepharid_ssu_PRJNA379597.fna",
        "analysis/markers/ssu_phylogeny/pr2.cdhit.reps"
    output:
        "analysis/markers/ssu_phylogeny/ssu.mafft"
    conda:
        "envs/mafft.yaml"
    threads:
        20
    shell:
        "cat {input} | mafft --thread {threads} --reorder --adjustdirection --auto - > {output}"

rule ssu_trim:
    input:
        "analysis/markers/ssu_phylogeny/ssu.mafft"
    output:
        "analysis/markers/ssu_phylogeny/ssu.trim"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -in {input} -out {output} -automated1"

rule collect_sintax:
    input:
        tsv = expand("analysis/markers/rna/ssu/{dataset}_{basename}.sintax.tsv", zip, dataset = datasets, basename = datasets.values()),
        datasets = "metadata/datasets.xlsx",
        clades = "output/clade_dist.csv"
    output:
        "output/sintax.svg"
    params:
        min_datasets = 25,
        min_rnas = 75,
        cutoff = 0.8,
        ranks = [ "f", "g", "s" ],
        include = "Eukaryota;",
        exclude = "Opisthokonta-Metazoa;",
        ingroup = "mgCCR1"
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/collect_sintax.R"
