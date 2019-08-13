#!/bin/bash

sed 's/chr//g' ucsc.known.genes.cds.hg19.bed | grep -v Y | grep -v X | grep -v hap | grep -v random | grep -v 17_ctg | grep -v 4_ctg \
  | grep -v Un_gl | grep -v 6_ssto | grep -v 19_gl | bedtools sort -i stdin -g human_g1k_v37_nochr.fasta.fai > ucsc.known.genes.cds.filtered.hg19.bed
