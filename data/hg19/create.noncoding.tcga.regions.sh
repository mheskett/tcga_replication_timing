#!/bin/bash

bedtools complement -i ucsc.known.genes.cds.filtered.hg19.bed -g human_g1k_v37_nochr.fasta.fai | \
  bedtools intersect -a stdin -b TCGA.segtabs.final.merged.sorted.coverage.9800.bed | 
 bedtools sort -i stdin -g human_g1k_v37_nochr.fasta.fai |  bedtools merge -d 100 -i stdin > noncoding.tcga.covered.regions.bed
