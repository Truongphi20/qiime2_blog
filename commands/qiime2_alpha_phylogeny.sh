python -m pdb -m q2cli phylogeny align-to-tree-mafft-fasttree \
  --i-sequences test_data/rep-seqs.qza \
  --o-alignment test_data/aligned-rep-seqs.qza \
  --o-masked-alignment test_data/masked-aligned-rep-seqs.qza \
  --o-tree test_data/unrooted-tree.qza \
  --o-rooted-tree test_data/rooted-tree.qza