# Deblur

## Introduction

Similar to [](./dada2.md), Deblur is a denoising method designed to correct errors occurring during sequencing and PCR amplification cycles, which otherwise limit the ability to perform fine-scale classification in 16S rRNA amplicon sequencing [@amirDeblurRapidlyResolves2017].

Instead of using traditional Operational Taxonomic Units (OTUs) or Amplicon Sequence Variants (ASVs) of DADA2, Deblur introduced a novel approach termed Sub Operational Taxonomic Units (sOTUs). 

There are two main steps in Deblur pipeline according to two commands mentioned in the tutorial.  

The command for preparation by quality filtering:
```
qiime quality-filter q-score \
  --i-demux demux.qza \
  --o-filtered-sequences demux-filtered.qza \
  --o-filter-stats demux-filter-stats.qza
```

And the core deblur command:
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

```{image} static/deblur_overview.png
:alt: dada_workflow
:height: 700px
:align: center
```

### Quality filtering

### Build index

### Deblur core

```{image} static/deblur_core.png
:alt: dada_workflow
:height: 500px
:align: center
```

### Create biom table

## Summary