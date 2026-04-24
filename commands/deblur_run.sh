python -m pdb -m q2cli denoise-16S \
  --i-demultiplexed-seqs /workspaces/qiime2_blog/test_data/demux-filtered.qza \
  --p-trim-length 120 \
  --p-sample-stats \
  --o-representative-sequences /workspaces/qiime2_blog/test_data/rep-seqs-deblur.qza \
  --o-table /workspaces/qiime2_blog/test_data/table-deblur.qza \
  --o-stats /workspaces/qiime2_blog/test_data/deblur-stats.qza