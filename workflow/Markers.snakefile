
from pandas import read_excel, read_csv
import warnings

# all datasets as { dataset: basename } dict
datasets = dict(zip(*glob_wildcards("datasets/{dataset}/{basename}.fna")))
markers ,= glob_wildcards("databases/eukaryota_odb10/hmms/{marker}.hmm")

# dataset metadata
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    metadata = read_excel("metadata/datasets.xlsx").set_index('analysis').T.to_dict('dict')

# dict of metatranscriptomic datasets
metat_datasets = { ds: datasets[ds] for ds in datasets.keys() if metadata[ds]['type'] == 'Metatranscriptome' }

rule all:
    input:
        expand("analysis/markers/{marker}.proteinortho.tsv", marker = markers)

rule cat_hmm:
    input:
        expand("databases/eukaryota_odb10/hmms/{marker}.hmm", marker = markers)
    output:
        "analysis/busco/markers.hmm"
    shell:
        "cat {input} > {output}"

rule getorf:
    input:
        "datasets/{dataset}/{basename}.fna"
    output:
        "analysis/getorf/{dataset}/{basename}.faa"
    params:
        minsize = 300 # min ORF size
    conda:
        "envs/emboss.yaml"
    shell:
        "getorf -minsize {params.minsize} -filter -find 0 {input} > {output}"

rule hmmsearch:
    input:
        fasta = "analysis/getorf/{dataset}/{basename}.faa",
        hmm = "analysis/busco/markers.hmm"
    output:
        txt = "analysis/hmmsearch/{dataset}/{basename}.txt",
        tblout = "analysis/hmmsearch/{dataset}/{basename}.tblout"
    conda:
        "envs/hmmer.yaml"
    threads:
        2
    shell:
        "hmmsearch --cpu {threads} --noali -o {output.txt} --tblout {output.tblout} {input.hmm} {input.fasta}"

rule hmmsearch_extract:
    input:
        faa = "analysis/getorf/{dataset}/{basename}.faa",
        fna = "datasets/{dataset}/{basename}.fna",
        tblout = "analysis/hmmsearch/{dataset}/{basename}.tblout",
        cutoffs = "databases/eukaryota_odb10/scores_cutoff"
    output:
        faa = expand("analysis/markers/{marker}/{{dataset}}_{{basename}}.faa", marker = markers),
        fna = expand("analysis/markers/{marker}/{{dataset}}_{{basename}}.fna", marker = markers)
    wildcard_constraints:
        dataset = '[^_]+'
    params:
        markers = markers
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/hmmsearch_extract.py"

rule proteinortho:
    input:
        expand("analysis/markers/{{marker}}/{dataset}_{basename}.fna", zip, dataset = metat_datasets, basename = metat_datasets.values())
    output:
        "analysis/markers/{marker}.proteinortho.tsv"
    params:
        ident = 90,
        evalue = 1e-20,
        prog = "blastn"
    shadow:
        "minimal"
    conda:
        "envs/proteinortho.yaml"
    shell:
        "proteinortho -project=myproject -cpus={threads} -identity={params.ident} -sim=0 -e {params.evalue} -p {params.prog} {input} && mv myproject.proteinortho.tsv {output}"



rule dload_busco:
    output:
        "analysis/busco/eukaryota"
    params:
        url = "https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz"
