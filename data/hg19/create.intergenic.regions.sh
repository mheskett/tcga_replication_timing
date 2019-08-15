#!/bin/bash

## want a file that is everything but introns and exons.

cat ucsc.introns.filtered.hg19.bed ucsc.exons.filtered.hg19.bed | \
  bedtools sort -i stdin -g human_g1k_v37_nochr.fasta.fai | \
  bedtools complement -i stdin -g human_g1k_v37_nochr.fasta.fai | \
  bedtools intersect -a stdin -b TCGA.segtabs.final.merged.sorted.coverage.9800.bed | \
  bedtools sort -i stdin -g human_g1k_v37_nochr.fasta.fai | \
  bedtools merge -i stdin > ucsc.not.i.e.tcga.covered.bed



