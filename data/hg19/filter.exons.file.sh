#!/bin/bash

sed 's/chr//g' ucsc.exons.hg19.bed | grep -v M | grep -v Y | grep -v X | grep -v hap | grep -v random | grep -v 17_ctg | grep -v 4_ctg \
  | grep -v Un_gl | grep -v 6_ssto | grep -v 19_gl | bedtools sort -i stdin -g human_g1k_v37_nochr.fasta.fai > ucsc.exons.filtered.hg19.bed
