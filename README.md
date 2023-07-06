## <a name="started"></a>Getting Started
```sh
# Align proteins to each genome
miniprot --outs=0.97 --no-cs -Iut16 genome1.fna proteins.faa > genome1.paf
miniprot --outs=0.97 --no-cs -Iut16 genome2.fna proteins.faa > genome2.paf

# Construct a pangene graph
pangene -a2 genome1.paf genome2.paf > graph.gfa

# Extract a subgraph around several genes
gfatools view -wl C4A,C4B -r3 graph.gfa > subgraph.gfa
# Visualize small subgraphs at https://lh3.github.io/gfatools/

# Check manpage
man ./pangene.1
```

## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [Usage](#usage)
  - [Preparing a protein set](#prep-aa)
  - [Aligning proteins to genomes](#align-aa)
  - [Constructing a pangene graph](#build-graph)
  - [Exploring a pangene graph](#explore-graph)
- [Limitations](#limit)

## <a name="intro"></a>Introduction

Pangene is a command-line tool to construct a pangenome gene graph. In this
graph, a node repsents a marker gene and an edge between two genes indicates
their genomic adjaceny on input genomes. Pangene takes the [miniprot][mp]
alignment between a protein set and multiple genomes and produces a graph in 
the GFA format. It attempts to reduce the redundancy in the input proteins and
filter spurious alignments while preserving close but non-identical paralogs.
The output graph can be visualized in generic GFA viewers such as
[BandageNG][bandage]. Users can also extract small subgraphs with
[gfatools][gfatools] and display with a simple [online GFA viewer][gfaview].
Prebuilt pangene graphs can be found at [DOI:10.5281/zenodo.8118576][zenodo].

Bacterial pangenome tools such as [panaroo][panaroo] often leverage gene graphs
to build bacterial pangenomes. Pangene is different in that it uses miniprot to
infer gene models. This makes pangene applicable to large Eukaryotic pangenomes
and robust to imperfect gene annotations.

Pangene is a **work-in-progress**. It may not be properly handling corner
cases during graph construction. Please create an issue if you see bugs or
questionable subgraphs.

## <a name="usage"></a>Usage

Pangene takes a list of protein-to-genome alignment as input. To generate
these alignments, you need to align the same set of proteins to multiple
genomes. How to choose the protein set can be tricky.

### <a name="prep-aa"></a>Preparing a protein set

For constructing a human pangene graph, the simplest choice is to use annotated
genes on GRCh38. It is highly recommended to name a protein sequence like
`RGPD6:ENSP00000512633.1` where `RGPD6` is the gene name and
`ENSP00000512633.1` is the protein identifier. Different isoforms of the same
gene can be present in the protein set. Pangene is designed to work with them.
In the output GFA, nodes are named after genes. You would want to use
human-readable gene names for visualization later.

Due to structural variations, some individuals may have genes distinct from the
gene annotations on the reference genome. In principle, it is preferred to
include structurally variable genes in the protein set. Nonetheless, such genes
are rare in the human genome. You can still get decent pangene graphs with
reference gene annotations only.

For a new species without good gene annotations, you may use protein annotations
from a closely related species. You may pool proteins from multiple closely
related species as well. Pangene aims to work with such input but this use case
has not been carefully evaluated.

Given a bacteria pangenome, you may predict genes with existing tools, cluster
them with CD-HIT or MMseqs2 and feed the representative protein of each cluster
to pangene. This apparently works for ~150 complete *Mycobacterium
tuberculosis* genomes but again, more evaluation is needed.

### <a name="align-aa"></a>Aligning proteins to genomes

Pangene currently only works with miniprot's PAF output. You may align proteins
to each genome with:
```sh
miniprot --outs=0.97 --no-cs -Iut16 genomeX.fna proteins.faa > genomeX.paf
```
For aligning proteins to bacterial genomes without splicing, add `-S` to the
command line above.

### <a name="build-graph"></a>Constructing a pangene graph

The following command-line constructs a pangene graph
```sh
pangene *.paf > graph.gfa
```
If the output graph is cluttered in the Bandage viewer, you may add option
`-a2` to filter out edges supported by a single genome. By default, pangene
filters out genes occurring in less than 5% of the genomes after deredundancy.
If you want to retain low-frequency genes, add `-p0` to disable the filter.

### <a name="explore-graph"></a>Exploring a pangene graph

You can visualize the entire pangene graph with the Bandage viewer. If you know
or find genes of interest, you can extract a subgraph with
```sh
gfatools view -wl C4A,C4B -r3 graph.gfa > subgraph.gfa
```
Here, `-w` tries to flip gene paths to the same orientation, `-l` specifies the
list of gene names and `-r` extracts nearby genes. If you put gene names in a
file `list.txt`, you may use
```sh
gfatools view -wl @list.txt -r3 graph.gfa > subgraph.gfa
```
You may visualize small subgraphs with the [online gfatools viewer][gfaview].
This viewer shows gene paths and counts their frequencies.

## <a name="limit"></a>Limitations

* In general, more testing needed.

* Pangene only works with [miniprot][mp]'s PAF output.

* In the output graph, arcs on W-lines may be absent from L-lines.

[mp]: https://github.com/lh3/miniprot
[bandage]: https://github.com/asl/BandageNG
[gfatools]: https://github.com/lh3/gfatools
[gfaview]: https://lh3.github.io/gfatools/
[panaroo]: https://github.com/gtonkinhill/panaroo
[asub]: https://github.com/lh3/asub
[zenodo]: https://zenodo.org/record/8118577
