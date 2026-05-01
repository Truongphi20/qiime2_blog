# Alpha diversity

## Introduction

Command for creating phylogenetic tree:
```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

Generating core matrices:
```bash
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1103 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir diversity-core-metrics-phylogenetic
```

Computing alpha diversity according a matrix:
```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity-core-metrics-phylogenetic/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization faith-pd-group-significance.qzv
```


## Workflow


## Summary