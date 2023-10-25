##################################
### Mycobacterium tuberculosis ###
##################################

# prodigal v2.6.3
# prokka v1.14.6
# panaroo v1.3.4
# miniprot v0.12

### gene finding ###

# training prodigal using the reference genome
prodigal -t prodigal_training_file -c -m -g 11 -p single -q -i H37Rv.fna > /dev/null

# prokka on all assemblies
ls *.fasta | sed s,.fasta,, | xargs -i echo prokka --prodigaltf prodigal_training_file --outdir {}.prokka {}.fasta \> {}.prokka.log 2\>\&1 | parallel -j7

# link gff files (as prokka name all of them as PROKKA_10192023.*)
ls *.LR.Asm.fasta|sed s,.LR.Asm.fasta,,|xargs -i echo ln -s {}.LR.Asm.prokka/PROKKA_10192023.gff {}.gff | sh

### pangene ###

# generate the protein set
cat H37Rv.faa *.prokka/*.faa | pigz -p8 > Mtb.all.faa.gz
cd-hit -i Mtb.all.faa.gz -o cluster -c .98 -M 0 -T 16 > cluster.log 2>&1
cat H37Rv.faa cluster | pigz -p8 > Mtb.faa.gz

# run miniprot
ls *.fasta | sed s,.fasta,, | xargs -i echo miniprot --outs=0.97 --no-cs -Iut8 {}.fasta ../Mtb.faa.gz \| gzip \> {}.paf.gz #| asub -j run-mp -M8 -n8

# generate pangene graphs
seqtk comp H37Rv.faa | cut -f1 > Mtb.H37Rv.txt  # prefer genes on the reference genome
pangene -P Mtb.H37Rv.txt -p.001 aln/*.paf.gz 2> Mtb-merge-p0.r177.gfa.log | gzip > Mtb-merge-p0.r177.gfa.gz &

### panaroo ###

# run panaroo
runlog panaroo -i *.gff -o panaroo.strict -t 16 --clean-mode strict > panaroo.strict.log 2>&1
runlog panaroo -i *.gff -o panaroo.strict-merge -t 16 --clean-mode strict --merge_paralogs > panaroo.strict-merge.log 2>&1

#############
### Human ###
#############

ls *.fa.gz | sed s,.fa.gz,, | xargs -i echo miniprot --outs=0.97 --no-cs -Iut16 {}.fa.gz ../HUMAN.faa.gz \| gzip \> {}.paf.gz|asub -j run-mp -M36 -n16
../pangene -a2 -I@gape.huamn.incl -X@gape.human.excl aln-human/[14]*.paf.gz 2> human-human98-r182.gfa.log | gzip > human-human98-r182.gfa.gz

##################
### Great apes ###
##################

# generate proteins
cat HUMAN.faa.gz PANTR.faa.gz PANPA.faa.gz GROGO.faa.gz PONAB.faa.gz > merged.faa.gz
ls *.fa.gz|sed s,.fa.gz,,|xargs -i echo miniprot --outs=0.97 --no-cs -Iut16 {}.fa.gz merged.faa.gz \| gzip \> {}.paf.gz | asub -j run-mp -M40 -n16
../pangene 1100_GRCh38.0.paf.gz 7031_chimpanzee.paf.gz 7032_bonobo.paf.gz 7033_gorilla.paf.gz 7034_orangutan.paf.gz > 5a.gfa
(zcat ../proteins/HUMAN.faa.gz; seqtk comp merged.faa.gz|cut -f1|grep -wFf <(../pangene.js gfa2matrix 5a.gfa|grep _|awk '$2==0'|cut -f1) | seqtk subseq -l80 merged.faa.gz -) | gzip > human+other.faa.gz
