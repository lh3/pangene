Release 1.1-r231 (22 February 2024)
-----------------------------------

Notable changes to pangene:

 * Improvement: more aggressive pseudogene filtering. This helps the great ape
   graph.

 * New feature: added option -E to filter out genes aligned with a single exon
   in all input genomes. Graphs constructed from multi-exon genes are usually
   cleaner.

Notable changes to pangene.js:

 * Improvement (call): a new bubble finding algorithm via local graph
   traversal. Added a filter to remove falsely idenrified bubbles.

 * New feature (call): output names of assemblies that support each allele.

 * New feature: added the getaa subcommand to parse Ensembl and GenCode
   annotations.

(22 February 2024, r231)
