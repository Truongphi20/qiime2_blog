# Deblur

## Introduction

Similar to [](./dada2.md), Deblur is a denoising method designed to correct errors occurring during sequencing and PCR amplification cycles, which otherwise limit the ability to perform fine-scale classification in 16S rRNA amplicon sequencing [@amirDeblurRapidlyResolves2017].

The command for quality filtering by quality scores:
```
qiime quality-filter q-score \
  --i-demux demux.qza \
  --o-filtered-sequences demux-filtered.qza \
  --o-filter-stats demux-filter-stats.qza
```


The deblur command:
```
qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 120 \
  --p-sample-stats \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --o-stats deblur-stats.qza
```


## Workflow

## Summary