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

The quality filtering process is performed independently for each sample and consists of three main steps: (1) scanning each read for a low-quality window, (2) truncating the FASTQ record based on the position of that window, and (3) tracing and labeling each record based on the results of the filter.

A low-quality window is identified as the first instance of consecutive bases with Phred quality scores below a specific threshold (the default is 4). Once this window is found, the read is truncated at that position. Each read is then assigned a status label to track its filtering outcome:

| Label   |   Meaning   |
| :----   | :---------- |
| **untruncated** |  The read was not truncated; the quality remained above the threshold throughout; do not contains ambiguous bases (`N`)        |
| **truncated**   |  The truncated read has accepted truncated fraction, and do not obtain ambiguous bases      |
| **short**       |  The read has truncated fraction (`truncated length` / `read length`) is greater than $0.75$         |
| **ambiguous**   |  The read without being truncated obtains ambiguous bases         |
| **truncated ambiguous** | The truncated read obtains ambiguous bases         |

(deblur-core)=
### Deblur core

```{image} static/deblur_core.png
:alt: dada_workflow
:height: 500px
:align: center
```

Qualify sequences **are trimed** to be equal in length (120 bp), sequences have length being shorter than the trim length are discard. **Deprelication** is performed by counting abundance of unique sequences and removing singletons using VSEARCH (v2.22.1) [@rognesVSEARCHVersatileOpen2016].

Following standardization, sequences were subjected to **multiple sequence alignment (MSA)** using MAFFT [@katohMAFFTMultipleSequence2013]. To accommodate the large-scale nature of the dataset, the PartTree algorithm was employed; this approach reduces the $O(N^2)$ complexity of pairwise comparisons to $O(N \log N)$ by recursively partitioning sequences based on their similarity to a subset of 'seed' sequences [@katohPartTreeAlgorithmBuild2007].

Metaphorically, **Launch Deblur** initializes a competition where each sequence is a competitor, and their "mana" represents their observed abundance. Starting with the most abundant, sequences duel one by one until the end, the survivors remain as sOTUs. In each match, the "damage points" dealt by a sequence to another represent the statistical number of error copies expected if the opponent were merely its own sequencing artifact (empirical algorithmic code can be found in [deblur.py](https://github.com/Truongphi20/qiime2_blog/blob/main/algorithm_reference/deblur.py)).

The **Chimera removal** step is performed by VSEARCH using the UCHIME *de novo* algorithm [@edgarUCHIMEImprovesSensitivity2011]. This approach operates on the assumption that 'parent' sequences coexist in the same FASTQ file as their chimeric artifacts. Any sequence with an abundance exceeding a specific threshold (default is 2) is considered a potential parent and stored in a local reference set. The algorithm processes sequences in order of decreasing abundance, if a query sequence is found to be a significant match-constructed from a combination of two parents in reference, it is flagged as a chimera and discarded.

### Create biom table

The denoised output from the [](#deblur-core) process across all samples is compiled into a single BIOM table. This table consists of a matrix (dimensions: $sOTUs \times samples$) where the cell values represent the remaining frequency (abundance) of each sequence. Note that samples containing zero reads are excluded from the matrix. Additionally, any sOTUs with a total cross-sample abundance falling below a specified threshold (default is 10) are discarded to filter out rare artifacts. Both the final BIOM table and the corresponding sOTU sequences are preserved for downstream analysis.

### Remove artifacts

SortMeRNA v2.0 [@kopylovaSortMeRNAFastAccurate2012] is employed to remove artifact sequences, which by default are composed of PhiX and sequencing adapters , see [artifact.fa](https://github.com/Truongphi20/qiime2_blog/blob/main/support_data/artifacts.fa). Artifact sequences that match the references are stored separately from the sOTU sequences. Any samples that become empty after this filtration are discarded, and the BIOM table and sequence file are updated accordingly.

## Summary