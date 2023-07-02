## Getting Started
```sh
# Align proteins to each genome
miniprot --outs=0.97 --no-cs -Iut16 genome1.fna proteins.faa > genome1.paf
miniprot --outs=0.97 --no-cs -Iut16 genome2.fna proteins.faa > genome2.paf

# Construct a pangene graph
pangene -a2 -e.5 genome1.paf genome2.paf > graph.gfa

# Extract a subgraph around several genes
gfatools view -wl C4A,C4B -r3 graph.gfa > subgraph.gfa

# Check manpage
man ./pangene.1
```
