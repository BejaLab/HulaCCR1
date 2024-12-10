HulaCCR1
========

This repository contains additional datasets and analysis workflow for the paper [HulaCCR1, a pump-like cation channelrhodopsin discovered in a lake microbiome](https://doi.org/10.1016/j.jmb.2024.168844)).

Files in the directory `metadata` (the input to the workflow):

* `BCCRs.faa` - representative bacteriorhodopsin-like cation channelrhodopsins (BCCRs)
* `datasets.xslx` - metadata for the datasets in which rhodopsins from the mgCCR1-3 clade (including HulaCCR1) were found
* `flanks.fna` - regions flanking genes coding for rhodopsin from the mgCCR1-3 clade
* `ingroup.faa` - all protein sequences of the collected rhodopsins from the mgCCR1-3 clade, includes fragments and redundant sequences
* `selected_ingroup.faa` - non-redundant set of representative proteins from the mgCCR1-3 clade

Files in the directory `output`:

* `clade_dist.csv` - distribution of the BCCRs in the different datasets containing ChRs from the mgCCR1-3 clade. CSV file with three columns: JGI GOLD analysis accession (or `HulaCCR1` for our Hula Lake dataset), BCCR (sub)clade, number of genes
* `pplacer.newick` - phylogenetic tree of the mgCCR1-3 clade with fragments placed on a [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) tree of representative sequences using [pplacer](https://matsen.fhcrc.org/pplacer/)
* `sintax.svg` - visualization of the distribution of the different taxa across the datasets containing genes from mgCCR1-3 clade based on 18S rRNA analysis with [sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html)
* `transposons_pplacer.newick` - phylogenetic tree of transposon sequences found in reference genomes and in the flanks of the contigs containing genes from the mgCCR1-3 clade. Fragments were placed on the tree with pplacer
* `transposons_subclades.csv` - summary of the transposon subclades identified in contigs with the genes from the mgCCR1-3 clade. For each representative sequence, numbers of fragments from each mgCCR1-3 subclade that cluster together with the representative are listed

Files in the directory `analysis`:

* `analysis/markers/rna/ssu/*.sintax.tsv` - raw results of the sintax analysis

The workflow itself is written in [Snakemake](https://snakemake.readthedocs.io/) and the dependencies are taken care of with [conda](https://anaconda.org/anaconda/conda). The code and the environment definitions can be found in the `workflow` directory.
