python -m pdb -m q2cli demux emp-single \
  --i-seqs test_data/emp-single-end-sequences.qza \
  --m-barcodes-file test_data/sample-metadata.tsv \
  --m-barcodes-column barcode-sequence \
  --o-per-sample-sequences test_data/demux.qza \
  --o-error-correction-details test_data/demux-details.qza