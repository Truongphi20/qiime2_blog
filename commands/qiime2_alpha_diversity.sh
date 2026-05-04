#!/bin/bash

OUTPUT="test_data/diversity-core-metrics-phylogenetic"
if [ -d $OUTPUT ]; then
  rm -r $OUTPUT
fi


python -m pdb -m q2cli diversity core-metrics-phylogenetic \
  --i-phylogeny test_data/rooted-tree.qza \
  --i-table test_data/table.qza \
  --p-sampling-depth 1103 \
  --m-metadata-file test_data/sample-metadata.tsv \
  --output-dir $OUTPUT


