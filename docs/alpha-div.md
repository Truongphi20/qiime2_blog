# Alpha diversity

## Introduction

Microbes has enormous effects on human life such as environment, heath and disease, and ecosystem. Acknowledging diversity of compositions in ecological units is essential fundamental input to answer ecological questions [@cassolKeyFeaturesGuidelines2025].

Depending on scale of ecological units, three levels of diversity measurement are: (1) Alpha diversity - describes the richness within a functional community. (2) Beta diversity - analyze about the differenciate between communities, (3) Gamma diversity - study on upper scale of diversity across geographic areas [@andermannEstimatingAlphaBeta2022].

In the QIIME 2 tutorial, focusing on alpha and beta diversity. Following steps are building phylogenetic tree, generating core matrices to analyse alpha and beta diversity.

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