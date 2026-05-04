# Diversity analysis

## Introduction

Microbes have a profound impact on human life, influencing the environment, health, disease, and ecosystems. Understanding the diversity and composition of these ecological units is fundamental to answering broader ecological questions [@cassolKeyFeaturesGuidelines2025].

Depending on the scale of the ecological units, diversity is measured at three levels [@andermannEstimatingAlphaBeta2022].

  - **Alpha diversity**: Describes the richness and evenness within a single community.
  - **Beta diversity**: Analyzes the differentiation between distinct communities.
  - **Gamma diversity**: Studies diversity at a larger landscape or geographic scale.

In the QIIME 2 tutorial, the focus is on alpha and beta diversity. The subsequent steps involve constructing a phylogenetic tree and generating core metrics matrices to facilitate these analyses.

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

## Workflow


## Summary