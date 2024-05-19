import os
from os import path
from pandas import read_excel, read_csv
import warnings

# all datasets as { dataset: basename } dict
datasets = dict(zip(*glob_wildcards("datasets/{dataset}/{basename}.fna")))

markers = {}
markers['busco'] ,= glob_wildcards("databases/eukaryota_odb10/hmms/{marker}.hmm")

markers['eukcc'] = []
eukcc_file = "databases/eukcc2_db_ver_1.2/db_protozoa/backbone/placement/protozoa_search.hmm"
with open(eukcc_file) as file:
    for line in file:
        if line.startswith("NAME"):
            markers['eukcc'].append(line.split()[1])

markers['rna'] = { "SSU": { "rnacentral": "URS0000726FAB", "rfam": "RF01960" }, "LSU": { "rnacentral": "URS0000ABD7EF", "rfam": "RF02543" } }
# dataset metadata
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    metadata = read_excel("metadata/datasets.xlsx").set_index('analysis').T.to_dict('dict')

# dict of metatranscriptomic datasets
metat_datasets = { ds: datasets[ds] for ds in datasets.keys() if metadata[ds]['type'] == 'Metatranscriptome' }

assemblages = { 'Freshwater': [ 'Freshwater_A', 'Freshwater_B', 'Freshwater_C', 'Freshwater_D' ], 'Marine': [ 'Marine_A' ] }
assemblages = { 'Freshwater': [ 'Freshwater_B' ], 'Marine': [ 'Marine_A' ] }

assemblages = { 'FreshwaterA': [ 'Freshwater_A' ], 'FreshwaterB': [ 'Freshwater_B' ], 'FreshwaterC': [ 'Freshwater_C' ], 'FreshwaterD': [ 'Freshwater_D' ], 'MarineA': [ 'Marine_A' ], 'MarineB': [ 'Marine_B' ] }

references = { 'Freshwater': 'HaHula1', 'Marine': 'Ga0314683' }
references = { 'FreshwaterA': 'HaHula1', 'FreshwaterB': 'HaHula1', 'FreshwaterC': 'HaHula1', 'FreshwaterD': 'HaHula1', 'MarineA': 'Ga0314683', 'MarineB': 'Ga0314683' }

db_datasets = {
    'rna': datasets,
    'busco': metat_datasets,
    'eukcc': metat_datasets
}

wildcard_constraints:
    dataset = "|".join(datasets)

rule all:
    input:
        # expand("analysis/markers/eukcc/{marker}_{assemblage}.selected.tsv", marker = markers['eukcc'], assemblage = assemblages),
        # expand("analysis/markers/busco/{marker}_{assemblage}.selected.tsv", marker = markers['busco'], assemblage = assemblages),
        expand("analysis/markers/{database}/{marker}_{assemblage}.proteinortho-selected.tsv", database = 'eukcc', marker = markers['eukcc'], assemblage = assemblages),
        expand("analysis/markers/{database}/{marker}_{assemblage}.proteinortho-selected.tsv", database = 'busco', marker = markers['busco'], assemblage = assemblages),
        expand("analysis/markers/{database}/{marker}_{assemblage}.proteinortho-selected.tsv", database = 'rna',   marker = markers['rna'],   assemblage = assemblages),

        expand("analysis/markers/{database}/{marker}_{assemblage}.uclust-selected.tsv", database = 'rna',   marker = markers['rna'],   assemblage = assemblages),
        expand("analysis/markers/{database}/{marker}_{assemblage}.uclust-selected.tsv", database = 'eukcc', marker = markers['eukcc'], assemblage = assemblages),
        expand("analysis/markers/{database}/{marker}_{assemblage}.uclust-selected.tsv", database = 'busco', marker = markers['busco'], assemblage = assemblages),

rule cat_busco_hmm:
    input:
        expand("databases/eukaryota_odb10/hmms/{marker}.hmm", marker = markers['busco'])
    output:
        "analysis/hmm/busco.hmm"
    shell:
        "cat {input} > {output}"

rule copy_busco_cutoffs:
    input:
        "databases/eukaryota_odb10/scores_cutoff"
    output:
        "analysis/hmm/busco_cutoffs.txt"

rule get_eukcc_cutoffs:
    input:
        "databases/eukcc2_db_ver_1.2/db_protozoa/backbone/placement/protozoa_search.hmm"
    output:
        "analysis/hmm/eukcc_cutoffs.txt"
    params:
        cut = 'GA'
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/hmm_cutoffs.py"

rule copy_eukcc_hmm:
    input:
        "databases/eukcc2_db_ver_1.2/db_protozoa/backbone/placement/protozoa_search.hmm"
    output:
        "analysis/hmm/eukcc.hmm"
    shell:
        "cp {input} {output}"

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
        hmm = "analysis/hmm/{database}.hmm"
    output:
        txt = "analysis/hmmsearch/{database}/{dataset}/{basename}.txt",
        tblout = "analysis/hmmsearch/{database}/{dataset}/{basename}.tblout"
    conda:
        "envs/hmmer.yaml"
    threads:
        2
    shell:
        "hmmsearch --cpu {threads} --noali -o {output.txt} --tblout {output.tblout} {input.hmm} {input.fasta}"

rule hmmsearch_eukprot:
    input:
        fasta = "databases/EukProt/proteins.fasta",
        hmm = "analysis/hmm/{database}.hmm"
    output:
        txt = "analysis/eukprot/{database}_hmmsearch.txt",
        tblout = "analysis/eukprot/{database}_hmmsearch.tblout"
    conda:
        "envs/hmmer.yaml"
    threads:
        20
    shell:
        "hmmsearch --cpu {threads} --noali -o {output.txt} --tblout {output.tblout} {input.hmm} {input.fasta}"

rule hmmsearch_eukprot_select:
    input:
        tblout = expand("analysis/eukprot/{database}_hmmsearch.tblout", database = [ 'eukcc', 'busco' ]),
        fasta = "databases/EukProt/proteins.fasta"
    output:
        "analysis/eukprot/markers.fasta"
    conda:
        "envs/kits.yaml"
    shell:
        "grep -hv '^#' {input.tblout} | cut -f1 -d' ' | sort -u | seqkit grep -f- -o {output} {input.fasta}"

rule get_unieuk:
    output:
        "analysis/eukprot/unieuk.tsv"
    params:
        url = "https://eukmap.unieuk.net/exports/unieuk/0.0.1-pre_release/unieuk-0.0.1-pre_release.tsv"
    shell:
        "wget -O {output} '{params.url}'"

rule get_eukprot_metadata:
    output:
        "analysis/eukprot/eukprot_metadata.tsv"
    params:
        url = "https://figshare.com/ndownloader/files/34436246" # EukProt_included_data_sets.v03.2021_11_22.txt
    shell:
        "wget -O {output} '{params.url}'"

rule eukprot_taxdb_taxonomy:
    input:
        fasta = "analysis/eukprot/markers.fasta",
        metadata = "analysis/eukprot/eukprot_metadata.tsv",
        unieuk = "analysis/eukprot/unieuk.tsv"
    output:
        names = "analysis/eukprot/names.dmp",
        nodes = "analysis/eukprot/nodes.dmp",
        mapping = "analysis/eukprot/mapping.txt"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/eukprot_taxdb.py"

rule eukprot_dummies:
    output:
        "analysis/eukprot/merged.dmp",
        "analysis/eukprot/delnodes.dmp"
    shell:
        "touch {output}"

rule eukprot_taxdb:
    input:
        fasta = "analysis/eukprot/markers.fasta",
        tax_dump = "databases/EukProt/taxdb",
        mapping = "analysis/eukprot/mapping.txt",
        nodes = "analysis/eukprot/nodes.dmp",
        names = "analysis/eukprot/names.dmp",
        merged = "analysis/eukprot/merged.dmp",
        delnodes = "analysis/eukprot/delnodes.dmp"
    output:
        "analysis/eukprot/markers_taxdb"
    conda:
        "envs/mmseqs.yaml"
    shell:
        "mmseqs createdb {input.fasta} {output} && mmseqs createtaxdb {output} tmp/ --ncbi-tax-dump $(dirname {input.nodes}) --tax-mapping-file {input.mapping}"

rule hmmsearch_extract_busco:
    input:
        faa = "analysis/getorf/{dataset}/{basename}.faa",
        fna = "datasets/{dataset}/{basename}.fna",
        tblout = "analysis/hmmsearch/busco/{dataset}/{basename}.tblout",
        cutoffs = "analysis/hmm/busco_cutoffs.txt"
    output:
        faa = expand("analysis/markers/busco/{marker}/{{dataset}}_{{basename}}.faa", marker = markers['busco']),
        fna = expand("analysis/markers/busco/{marker}/{{dataset}}_{{basename}}.fna", marker = markers['busco'])
    params:
        markers = markers['busco'],
        cutoff_factor = 0.5
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/hmmsearch_extract.py"

rule hmmsearch_extract_eukcc:
    input:
        faa = "analysis/getorf/{dataset}/{basename}.faa",
        fna = "datasets/{dataset}/{basename}.fna",
        tblout = "analysis/hmmsearch/eukcc/{dataset}/{basename}.tblout",
        cutoffs = "analysis/hmm/eukcc_cutoffs.txt"
    output:
        faa = expand("analysis/markers/eukcc/{marker}/{{dataset}}_{{basename}}.faa", marker = markers['eukcc']),
        fna = expand("analysis/markers/eukcc/{marker}/{{dataset}}_{{basename}}.fna", marker = markers['eukcc'])
    params:
        markers = markers['eukcc'],
        cutoff_factor = 0.5
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/hmmsearch_extract.py"

def non_empty_files(wildcards, input):
    return [ f for f in input if os.stat(f).st_size > 0 ]

rule proteinortho:
    input:
        expand("analysis/markers/{{database}}/{{marker}}/{dataset}_{basename}.fna", zip, dataset = metat_datasets, basename = metat_datasets.values())
    output:
        "analysis/markers/{database}/{marker}.proteinortho.tsv"
    wildcard_constraints:
        database = 'busco|eukcc'
    params:
        ident = 90,
        evalue = 1e-20,
        prog = "blastn",
        cov = 20,
        conn = 0.05,
        files = non_empty_files
    shadow:
        "minimal"
    conda:
        "envs/proteinortho.yaml"
    threads:
        10
    shell:
        "proteinortho -selfblast -project=p -cpus={threads} -identity={params.ident} -cov={params.cov} -conn={params.conn} -e={params.evalue} -p={params.prog} {params.files} && mv p.proteinortho.tsv {output}"

rule proteinortho_rna:
    input:
        expand("analysis/markers/rna/{{marker}}/{dataset}_{basename}.fna", zip, dataset = datasets, basename = datasets.values())
    output:
        "analysis/markers/rna/{marker}.proteinortho.tsv"
    params:
        ident = 95,
        evalue = 1e-20,
        prog = "blastn",
        cov = 20,
        conn = 0.05,
        files = non_empty_files
    shadow:
        "minimal"
    conda:
        "envs/proteinortho.yaml"
    threads:
        10
    shell:
        "proteinortho -selfblast -project=p -cpus={threads} -identity={params.ident} -cov={params.cov} -conn={params.conn} -e={params.evalue} -p={params.prog} {params.files} && mv p.proteinortho.tsv {output}"

rule uclust_input:
    input:
        lambda w: expand("analysis/markers/{{database}}/{{marker}}/{dataset}_{basename}.fna", zip, dataset = db_datasets[w.database], basename = db_datasets[w.database].values())
    output:
        "analysis/markers/{database}/{marker}_all.fasta"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit sort -rl {input} -o {output}"

rule uclust:
    input:
        "analysis/markers/{database}/{marker}_all.fasta"
    output:
        centroids = "analysis/markers/{database}/{marker}.uclust.fasta",
        clusters = "analysis/markers/{database}/{marker}.uclust.txt"
    # not in conda
    params:
        id = lambda w: 0.95 if w.database == 'rna' else '0.9'
    shell:
        "usearch -cluster_fast {input} -strand both -id {params.id} -centroids {output.centroids} -uc {output.clusters}"

rule uclust_wide:
    input:
        "analysis/markers/{database}/{marker}.uclust.txt"
    output:
        "analysis/markers/{database}/{marker}.uclust.tsv"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/uclust_wide.py"

rule filter_genes:
    input:
        tsv = "analysis/markers/{database}/{marker}.{strategy}.tsv",
        clades = "output/clade_dist.csv"
    output:
        "analysis/markers/{database}/{marker}_{assemblage}.{strategy}-selected.tsv"
    params:
        min_datasets = 0,
        assemblage = lambda w: assemblages[w.assemblage],
        reference = lambda w: references[w.assemblage],
        datasets = lambda w: db_datasets[w.database]
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/filter_genes.py"

rule make_mmseqs_db:
    input:
        "analysis/markers/{database}/{marker}/{filename}.faa"
    output:
        "analysis/markers/{database}/{marker}/{filename}_db"
    conda:
        "envs/mmseqs.yaml"
    shell:
        "mmseqs createdb {input} {output}"

rule taxonomy:
    input:
        query = "analysis/markers/{database}/{marker}/{filename}_db",
        taxdb = "analysis/eukprot/markers_taxdb"
    output:
        "analysis/markers/{database}/{marker}/{filename}_tax"
    params:
        maj = 0.5,
        mode = 2
    conda:
        "envs/mmseqs.yaml"
    threads:
        4
    shell:
        "mmseqs taxonomy {input.query} {input.taxdb} {output} tmp/ --lca-mode {params.mode} --majority {params.maj} --threads {threads}"

rule taxonomy_tsv:
    input:
        query = "analysis/markers/{database}/{marker}/{filename}_db",
        result = "analysis/markers/{database}/{marker}/{filename}_tax"
    output:
        "analysis/markers/{database}/{marker}/{filename}_tax.tsv"
    conda:
        "envs/mmseqs.yaml"
    shell:
        "mmseqs createtsv {input.query} {input.result} {output}"

rule collect_taxonomy:
    input:
        lambda w: expand(expand("analysis/markers/{database}/{marker}/{{dataset}}_{{basename}}_tax.tsv", marker = markers[w.database], database = w.database), zip, dataset = metat_datasets, basename = metat_datasets.values())
    output:
        "analysis/markers/{database}_taxonomy.tsv"
    conda:
        "envs/kits.yaml"
    shell:
        "awk 'NR==1||FNR>1' {input} > {output}"

rule dload_rnacentral:
    output:
        "analysis/rna/{rna}.fna"
    params:
        url = lambda w: "https://rnacentral.org/api/v1/rna/{accession}.fasta".format(accession = markers['rna'][w.rna]['rnacentral'])
    shell:
        "wget -O {output} {params.url}"

rule dload_rfam:
    output:
        expand("analysis/rna/{acc}.cm", acc = [ marker['rfam'] for marker in markers['rna'].values() ])
    params:
        url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz",
        dirname = lambda w, output: path.dirname(output[0]),
        files = lambda w, output: list(map(path.basename, output))
    shell:
        "wget -O- '{params.url}' | tar xfz - -C {params.dirname} {params.files}"

rule makeblastdb_rna:
    input:
        "analysis/rna/{rna}.fna"
    output:
        "analysis/rna/{rna}.fna.ndb"
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule blast_rna:
    input:
        db = "analysis/rna/{marker}.fna",
        ndb = "analysis/rna/{marker}.fna.ndb",
        query = "datasets/{dataset}/{basename}.fna"
    output:
        "analysis/blastn/rna/{marker}/{dataset}_{basename}.txt"
    params:
        evalue = 1e-10
    threads:
        4
    shell:
        "blastn -query {input.query} -db {input.db} -out {output} -outfmt 6 -evalue {params.evalue} -max_target_seqs 1000000000 -num_threads {threads}"

rule extract_rna:
    input:
        blast = "analysis/blastn/rna/{marker}/{dataset}_{basename}.txt",
        fna = "datasets/{dataset}/{basename}.fna"
    output:
        "analysis/markers/rna/{marker}/{dataset}_{basename}_full.fna"
    params:
        m = 150
    shell:
        "[ -s {input.blast} ] && cut -f1 {input.blast} | sort -u | seqkit grep -f- {input.fna} | seqkit seq -gm {params.m} | seqkit replace -p ^ -r '{wildcards.dataset}|' > {output} || touch {output}"

rule cmsearch:
    input:
        fasta = "analysis/markers/rna/{marker}/{dataset}_{basename}_full.fna",
        cm = lambda w: "analysis/rna/{acc}.cm".format(acc = markers['rna'][w.marker]['rfam'])
    output:
        "analysis/markers/rna/{marker}/{dataset}_{basename}.tblout"
    conda:
        "envs/infernal.yaml"
    threads:
        4
    shell:
        "[ -s {input.fasta} ] && cmsearch --cpu {threads} --noali --tblout {output} {input.cm} {input.fasta} || touch {output}"

rule trim_rna:
    input:
        tblout = "analysis/markers/rna/{marker}/{dataset}_{basename}.tblout",
        fasta = "analysis/markers/rna/{marker}/{dataset}_{basename}_full.fna"
    output:
        "analysis/markers/rna/{marker}/{dataset}_{basename}.fna"
    params:
        min_len = 150
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/trim_rna.py"
