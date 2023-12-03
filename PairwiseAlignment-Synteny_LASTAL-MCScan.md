## Replace sequence headers of HiC assemblies with the identified chromosomes for ManJav1 and ManPen1. Change header to numbers for PhaTri1 to make visualization and processing easier.
# create a tab-delimited key with original sequence name is the first column and the new name in the second column. Make sure to format assembly keys using awk:
# make sure to standarized the chromosome sequence names and the assembly name/code throughout this script to prevent any errors in process.

awk -v OFS="\t" '{print $1,$2}' assembly_key.txt > <assembly_name>_ChrName.txt

seqkit replace -w0 -p '(HiC_scaffold_[0-9]+)' -r '{kv}${2}' -k ManJav1_ChrName.txt ManJav1_RM/ManJav1.0_cscaffolds.fasta.masked | seqkit sort -n -N -w0 -o ManJav1.chromosomes.fasta -
seqkit replace -w0 -p '(HiC_scaffold_[0-9]+)' -r '{kv}${2}' -k ManPen1_ChrName.txt ManPen1_RM/ManPen1.0_cscaffolds.fasta.masked | seqkit sort -n -N -w0 -o ManPen1.chromosomes.fasta -
seqkit replace -w0 -p '(HiC_scaffold_[0-9]+)' -r '{kv}${2}' -k PhaTri1_ChrName.txt PhaTri1_RM/PhaTri1.0_cscaffolds.fasta.masked | seqkit sort -n -N -w0 -o PhaTri1.chromosomes.fasta -

## create LASTAL database for the reference that will be aligned to in all pairwise comparisons
# the lastal command that is uncommented is a faster version while the commented one is more accurate but takes longer to run
lastdb -P14 -uNEAR -cR11 manJav.db ManJav1.chromosomes.fasta

last-train -P14 --revsym --matsym --gapsym -E0.05 -C2 manJav.db PhaTri1.chromosomes.fasta > ManJav1.PhaTri1.mat
lastal -P14 -E0.05 --split-f=MAF+ -p ManJav1.PhaTri1.mat PhaTri1.chromosomes.fasta | maf-swap - > ManJav1.PhaTri1.many-to-one.maf
#lastal -P14 -m50 -E0.05 --split-f=MAF+ -p ManJav1.PhaTri1.mat PhaTri1.chromosomes.fasta | maf-swap - > ManJav1.PhaTri1.many-to-one.maf
last-split -r ManJav1.PhaTri1.many-to-one.maf | maf-swap - > ManJav1.PhaTri1.maf

last-train -P14 --revsym --matsym --gapsym -E0.05 -C2 manJav.db ManPen1.chromosomes.fasta > ManJav1.ManPen1.mat
lastal -P14 -E0.05 --split-f=MAF+ -p ManJav1.ManPen1.mat manJav.db ManPen1.chromosomes.fasta | maf-swap - > ManJav1.ManPen1.many-to-one.maf
#lastal -P14 -m50 -E0.05 --split-f=MAF+ -p ManJav1.ManPen1.mat manJav.db ManPen1.chromosomes.fasta | maf-swap - > ManJav1.ManPen1.many-to-one.maf
last-split -r ManJav1.ManPen1.many-to-one.maf | maf-swap - > ManJav1.ManPen1.maf

## Modify maf files and generate record files for MCScan (adapted from Yin et al. 2021 github repository)
perl maf.rename.species.S.pl ManJav1.ManPen1.maf ManJav1 ManPen1 ManJav1.ManPen1.swap.name.maf
perl maf.rename.species.S.pl ManJav1.PhaTri1.maf ManJav1 PhaTri1 ManJav1.PhaTri1.swap.name.maf
maf2region_lzs.pl ManJav1.ManPen1.swap.name.maf ManJav1.ManPen1
maf2region_lzs.pl ManJav1.PhaTri1.swap.name.maf ManJav1.PhaTri1

## These next commands generate all the needed files for MCScan from the record files above.
# Again, remember to standardize the sequence headers and assembly names/codes
# The alignment size filter for these commands is 14,000bp which you can file in the awk pipe at the end of the command. Modify this name to change the alignment size filter.
grep -v 'sv' ManJav1.ManPen1.record | sed '/^\s*$/d;/\*/d' | awk '$3-$2>=14000{print$0}' > ManJav1.ManPen1.filter.record
grep -v 'sv' ManJav1.PhaTri1.record | sed '/^\s*$/d;/\*/d' | awk '$3-$2>=14000{print$0}' > ManJav1.PhaTri1.filter.record

awk '{print $1"_"$2"\t"$1"_"$3"\t"$6"_"$7"\t"$6"_"$8"\t500\t"$9}' ManJav1.ManPen1.filter.record | sed 's/ManJav1.chr//g;s/ManPen1.chr//g' > ManJav1.ManPen1.simple   
awk '{print $1"\t"$2-1"\t"$2+1"\t"$1"_"$2"\t0\t"$4"\n"$1"\t"$3-1"\t"$3+1"\t"$1"_"$3"\t0\t"$4}' ManJav1.ManPen1.filter.record > ManJav1.ManPen1.ManJav1.bed  
awk '{print $6"\t"$7-1"\t"$7+1"\t"$6"_"$7"\t0\t"$9"\n"$6"\t"$8-1"\t"$8+1"\t"$6"_"$8"\t0\t"$9}' ManJav1.ManPen1.filter.record > ManJav1.ManPen1.ManPen1.bed

awk '{print $1"_"$2"\t"$1"_"$3"\t"$6"_"$7"\t"$6"_"$8"\t500\t"$9}' ManJav1.PhaTri1.filter.record | sed 's/ManJav1.chr//g;s/PhaTri1.chr//g' > ManJav1.PhaTri1.simple
awk '{print $1"\t"$2-1"\t"$2+1"\t"$1"_"$2"\t0\t"$4"\n"$1"\t"$3-1"\t"$3+1"\t"$1"_"$3"\t0\t"$4}' ManJav1.PhaTri1.filter.record > ManJav1.PhaTri1.ManJav1.bed 
awk '{print $6"\t"$7-1"\t"$7+1"\t"$6"_"$7"\t0\t"$9"\n"$6"\t"$8-1"\t"$8+1"\t"$6"_"$8"\t0\t"$9}' ManJav1.PhaTri1.filter.record > ManJav1.PhaTri1.PhaTri1.bed 

cat ManJav1.ManPen1.ManJav1.bed ManJav1.PhaTri1.ManJav1.bed | sed 's/ManJav1.chr//g' > ManJav1.bed
cat ManJav1.ManPen1.ManPen1.bed | sed 's/ManPen1.chr//g' > ManPen1.bed
cat ManJav1.PhaTri1.PhaTri1.bed | sed 's/PhaTri1.chr//g' > PhaTri1.bed

## examples of seqids file and layout files
# # pangolin_seqids
# 1,3,2,4,5,8,11,10,17,7,12,6,9,15,18,13,14,16,19,X
# 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X
# 29,20,1,33,35,53,23,31,21,5,8,55,10,14,3,17,50,54,9,37,57,7,44,26,40,12,34,6,56,24,41,11,27,13,46,22,48,45,42,51,4,25,36,49,32,18,19,39,38,15,43,16,30,52,47,28,X
# 
# # layout_pangolins
# # y, xstart, xend, rotation, color, label, va,  bed
#  .7,     .1,    .8,       0,   red, ManPen1, top, ManPen1.bed_final
#  .5,     .1,    .8,       0,  blue, ManJav1, top, ManJav1.bed
#  .3,     .1,    .8,	  0, green, PhaTri1, bottom, PhaTri1.bed_final
# # edges
# e, 0, 1, ManJav1.ManPen1.simple_highlightFinal
# e, 1, 2, ManJav1.PhaTri1.simple_highlightFinal

## run jcvi mcscan to generate an image plot that can be further edited as you need, you will have to run this again after you identify any chromosomes that need to be uninverted
python -m jcvi.graphics.karyotype pangolin_seqids layout_pangolins --notex --figsize=13x11 --dpi=600 --format=pdf -o Pangolin_karyotype.pdf

## example of the key file for the two shell loops below
# # phatri1.key
# original	20
# 20	33
# 33	35
# 35	31
# 31	55
# 55	17
# 17	50
# 50	9
# 9	7
# 7	6
# 6	24
# 24	11
# 11	13
# 13	51
# 51	4
# 4	36
# 36	49
# 49	32
# 32	18
# 18	39
# 39	15
# 15	43
# 43	30
# 30	final
# 
# # color.key
# 2	purple
# 3	teal
# X	olive

## for inverting chromosomes create a two-column text file with  "original" in the first row and then going down a list of the chromosomes that are being inverted. In the second column, start value in the first row and proceed downwards. At the end put "final".
while read x y \n
do echo ${x} \n 
python Assembly/MuntjacAnalysis/mcscan_invert_chr.py PhaTri1.bed_${x} phatri1.renamed.fa.fai ${y} > PhaTri1.bed_${y} \n
done < phatri1.key

## create then same formatted two-column text file as before
cp ManJav1.ManPen1.simple ManJav1.ManPen1.simple_highlight-${x}
read while x y z
do
echo "Highlighting ${x}"
sed "s/^${x}_/${y}*${x}_/g' ManJav1.ManPen1.simple_highlight-${x} > ManJav1.ManPen1.simple_highlight-${y}
done < color.key

# generate finalized synteny plot
python -m jcvi.graphics.karyotype pangolin_seqids layout_pangolins --notex --figsize=13x11 --dpi=600 --format=pdf -o Pangolin_karyotype.pdf
