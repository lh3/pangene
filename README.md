## Getting Started
```sh
# Align proteins to each genome
miniprot --outs=0.97 --no-cs -Iut16 genome1.fna proteins.faa > genome1.paf
miniprot --outs=0.97 --no-cs -Iut16 genome2.fna proteins.faa > genome2.paf

# Construct a pangene graph
pangene -a2 -e.5 genome1.paf genome2.paf > graph.gfa

# Extract a subgraph around several genes
gfatools view -wl C4A,C4B -r3 graph.gfa > subgraph.gfa
# Visualize small subgraphs at https://lh3.github.io/gfatools/

# Check manpage
man ./pangene.1
```

## Introduction

Pangene is a command-line tool to construct a pangenome gene graph. In this
graph, a node repsents a marker gene and an edge between two genes indicates
their genomic adjaceny on input genomes. Pangene takes the [miniprot][mp]
alignment between one protein set and multiple genomes and outputs a graph in
the GFA format. It attempts to reduce the redundancy in input proteins and to
filter spurious alignments while preserving close but non-identical paralogs.
The output graph can be visualized in generic GFA viewers such as
[BandageNG][bandage]. Users can also extract small subgraphs with
[gfatools][gfatools] and display with a simple [online GFA viewers][gfaview]
which have recently been updated to support the pangene output.

Bacterial pangenome tools such as [panaroo][panaroo] often leverage gene graphs
to build bacterial pangenomes. Pangene is different in that it uses miniprot to
infer gene models and works with large Eukaryotic pangenomes. Pangene has been
used to build a gene graph for ~150 complete *Mycobacterium tuberculosis*
genomes with CD-HIT clusters as input, but the quality of the output graph has
not been carefully evaluated.

Pangene is a **work-in-progress**. It may not be properly handling subtle
details during graph construction. Please create an issue if you see bugs or
questionable subgraphs.

## Limitations

* In general, more testing needed.

* Pangene only works with [miniprot][mp]'s PAF output.

* In the output graph, arcs on W-lines may be absent from L-lines.

[mp]: https://github.com/lh3/miniprot
[bandage]: https://github.com/asl/BandageNG
[gfatools]: https://github.com/lh3/gfatools
[gfaview]: https://lh3.github.io/gfatools/
[panaroo]: https://github.com/gtonkinhill/panaroo
