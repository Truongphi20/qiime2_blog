python -m pdb -m q2cli dada2 denoise-single \
  --i-demultiplexed-seqs test_data/demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 120 \
  --o-representative-sequences test_data/rep-seqs.qza \
  --o-table test_data/table.qza \
  --o-denoising-stats test_data/denoising-stats.qza \
  --o-base-transition-stats test_data/base-transition-stats.qza