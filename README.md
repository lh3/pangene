## <a name="started"></a>Getting Started
```sh
# Check prebuilt graphs at https://pangene.bioinweb.org

# Install pangene
git clone https://github.com/lh3/pangene
cd pangene && make

# Alternatively, download the precompiled binaries for arm-mac and x64-linux
curl -L https://github.com/lh3/pangene/releases/download/v1.1/pangene-1.1-bin.tar.bz2|tar jxf -

# The C4 example with provided alignment
./pangene test/C4/*.paf.gz > C4.gfa        # generate the graph
k8 pangene.js call C4.gfa > C4.bubble.txt  # identify bubbles

# Deploy the GFA server on the C4 example; require pangene-1.1-bin
cd pangene-1.1-bin                         # run gfa-server in this directory
bin_arm64-mac/gfa-server /path/to/C4.gfa   # or use bin_linux-x64 on x64-linux
# open "http://127.0.0.1:8000/view?gene=C4A,C4B" in a web browser

# Deploy the GFA server on the human graph
bin_arm64-mac/gfa-server -d html data/*.gfa.gz 2> server.log &
# open "http://127.0.0.1:8000" in a web browser

# Align proteins to each genome (general use cases; no example data)
miniprot --outs=0.97 --no-cs -Iut16 genome1.fna proteins.faa > genome1.paf
miniprot --outs=0.97 --no-cs -Iut16 genome2.fna proteins.faa > genome2.paf

# Construct a pangene graph
pangene genome1.paf genome2.paf > graph.gfa

# Check manpage
man ./pangene.1
```

## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [Graph Construction](#build)
  - [Preparing a protein set](#prep-aa)
  - [Aligning proteins to genomes](#align-aa)
  - [Constructing a pangene graph](#build-graph)
  - [Analyzing a graph](#analyze)
- [Graph Visualization](#visual)
- [Limitations](#limit)

## <a name="intro"></a>Introduction

Pangene is a command-line tool to construct a pangenome gene graph. In this
graph, a node repsents a marker gene and an edge between two genes indicates
their genomic adjaceny on input genomes. Pangene takes the [miniprot][mp]
alignment between a protein set and multiple genomes and produces a graph in 
the GFA format. It attempts to reduce the redundancy in the input proteins and
filter spurious alignments while preserving close but non-identical paralogs.
The output graph can be visualized in generic GFA viewers such as
[BandageNG][bandage] or via a [web interface](@visual). Users can explore local
human subgraphs at a [public server][server]. Prebuilt pangene graphs can be
found at [DOI:10.5281/zenodo.8118576][zenodo].

<!--
Bacterial pangenome tools such as [panaroo][panaroo] often leverage gene graphs
to build bacterial pangenomes. Pangene is different in that it uses miniprot to
infer gene models. This makes pangene applicable to large Eukaryotic pangenomes
and robust to imperfect gene annotations.
-->

## <a name="build"></a>Graph Construction

Pangene takes a list of protein-to-genome alignment as input. To generate
these alignments, you need to align the same set of proteins to multiple
genomes. How to choose the protein set can be tricky.

### <a name="prep-aa"></a>Preparing a protein set

For constructing a human pangene graph, the simplest choice is to use annotated
genes in GRCh38. It is highly recommended to name a protein sequence like
`RGPD6:ENSP00000512633.1` where `RGPD6` is the gene name and
`ENSP00000512633.1` is the unique protein identifier. In the output GFA, nodes
are named after genes, so you would want to use human-readable gene names for
visualization later. You may use the following command line to extract protein
sequences from Ensembl or GenCode annotation:
```sh
k8 pangene.js getaa gene-anno.gtf protein-seq.faa > proteins.faa
```

With pangene, different isoforms or diverged alleles of the same gene can be
present in the protein set, though in practice, we find selecting the canonical
isoform per gene tends to give a cleaner graph probably possibly due to
annotation errors among rare isoforms. For the GenCode annotation, use `getaa
-c` to extract canonical isoforms only.

For a new species without good gene annotation, you may use protein annotations
from a closely related species. You may pool proteins from multiple closely
related species as well. Pangene aims to work with such input but this use case
has not been thoroughly carefully evaluated. Given a bacteria pangenome of the
same species, you may predict genes with existing tools, cluster them with
CD-HIT or MMseqs2 and feed the representative protein in each cluster to
pangene.

### <a name="align-aa"></a>Aligning proteins to genomes

Pangene currently only works with miniprot's PAF output. We usually use the
following command line:
```sh
miniprot --outs=0.97 --no-cs -Iut16 genomeX.fna proteins.faa > genomeX.paf
```
For bacterial data, add `-S` to disable splicing.

### <a name="build-graph"></a>Constructing a pangene graph

The following command-line constructs a pangene graph
```sh
pangene *.paf > graph.gfa
```
If the output graph is cluttered in the Bandage viewer, you may add option
`-a2` to filter out edges supported by a single genome. By default, pangene
filters out genes occurring in less than 5% of the genomes after deredundancy.
If you want to retain low-frequency genes, add `-p0` to disable the filter.

### <a name="analyze"></a>Analyzing a graph

The GFA file is the master output. You can extract various information from
this file. You may find local gene-level variations with
```sh
k8 pangene.js call graph.gfa > bubble.txt
```
or get the presence/absence of each gene with
```sh
k8 pangene.js gfa2matrix graph.gfa > gene_presence_absence.Rtab
```

## <a name="visual"></a>Graph Visualization

You can look at the entire graph in the Bandage GFA viewer. If you are
interested in a particular gene, it is best to set up gfa-server which is part
of [gfatools][gfatools]. [Here][server] is a public server for human genes.
You can deploy this server on your machine with
```sh
curl -L https://github.com/lh3/pangene/releases/download/v1.1/pangene-1.1-bin.tar.bz2|tar jxf -
cd pangene-1.1-bin
bin_arm64-mac/gfa-server -d html data/*.gfa.gz 2> server.log # for Mac
```
Then you can open link `http://127.0.0.1:8000/` in your browser, type gene
names and visualize a local subgraph around the desired genes.

## <a name="limit"></a>Limitations

* Pangene only works with [miniprot][mp]'s PAF output.

* In the output graph, arcs on W-lines may be absent from L-lines.

[mp]: https://github.com/lh3/miniprot
[bandage]: https://github.com/asl/BandageNG
[gfatools]: https://github.com/lh3/gfatools
[gfaview]: https://lh3.github.io/gfatools/
[panaroo]: https://github.com/gtonkinhill/panaroo
[asub]: https://github.com/lh3/asub
[zenodo]: https://doi.org/10.5281/zenodo.8118576
[server]: https://pangene.bioinweb.org
