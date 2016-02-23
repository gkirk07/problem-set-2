#! /usr/bin/env bash

#Question 1: Use BEDtools intersect to identify the size of the largest
#overlap between CTCF and H2K4me3 locations.

datasets=/Users/gregkirkpatrick/devel/data-sets

tfbs=$datasets/bed/encode.tfbs.chr22.bed.gz
h3k4=$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz

gzcat $tfbs | awk '$4 == "CTCF"' > ctcf-peaks.bed

answer_1=$(bedtools intersect -wo -a ctcf-peaks.bed -b $h3k4 | awk '{print $NF}' | sort -nr \
| head -n1)

echo "answer-1: $answer_1"

#Question 2:Use BEDtools to calculate the GC content of nucleotides
#19,000,000 to 19,000,500 on chr22 of hg19 genome build. Report the GC
#content as a fraction (e.g., 0.50).

echo -e 'chr22\t19000000\t19000500' > region.bed
answer_2=$(bedtools nuc -fi $datasets/fasta/hg19.chr22.fa -bed region.bed | awk '{print $5}' | grep -v "^5")

echo "answer-2: $answer_2"


#Question 3: Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
#interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz.

hela=$datasets/bedtools/ctcf.hela.chr22.bg.gz

answer_3=$(bedtools map -a <(gzcat $tfbs | awk '$4 == "CTCF"') -b $hela -c 4 -o mean | sort -k5rn | head -n1 | awk '{print $3-$2}')

echo "answer-3: $answer_3"

#Question 4: Use BEDtools to identify the gene promoter (1000 bp upstream
# of a TSS) with the highest median signal in ctcf.hela.chr22.bg.gz.
# Report the gene name

tss=$datasets/bed/tss.hg19.chr22.bed.gz
bedgraph=$datasets/bedtools/ctcf.hela.chr22.bg.gz
genome=$datasets/bedtools/hg19.genome
answer_4=$(gzcat $tss | bedtools flank -s -l 1000 -r 0 -g $genome \
| bedtools sort -i - \
| bedtools map -a - -b $bedgraph -c 4 -o median \
| sort -k7rn \
| head -n1 | cut -f4)

echo "answer-4: $answer_4"

#Question 5: Use BEDtools to identify the longest interval on chr22 that
#is not covered by genes.hg19.bed.gz.  Report the interval (chr1:100-500)

genes=$datasets/bed/genes.hg19.bed.gz

answer_5=$(bedtools complement -i $genes -g $genome | grep -w "chr22" \
| awk '{print $1, $2, $3, $3 - $2}' \
| sort -k4rn \
| head -n1 \
| awk 'BEGIN {OFS=""} {print $1, ":", $2, "-", $3}')

echo "answer-5: $answer_5"

#Question 6: What's the name of the closest gene on chromosome 22 between
#nucleotides 19,000,000 and 19,000,500?

gzcat $genes | awk '$1 == "chr22"' | bedtools sort -i - > chr22genes.bed
answer_6=$(bedtools closest -a region.bed -b chr22genes.bed | cut -f7)

echo "answer-6: $answer_6"


