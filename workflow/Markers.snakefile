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
