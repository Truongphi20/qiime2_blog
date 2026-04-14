# Demultiplexing

## Introduction

"Multiplexing" is a sequencing technique, which pools many samples in the same sequecning batch, where samples are distinguished by barcode sequences. It helps enhance throughput, optimize cost, and simplifiy analysis [@SampleMultiplexingMultiplex]. In contrast, "demultiplexing" resolves that which sample that sequencing read belonging. 

The command in the tutorial:
```bash
qiime demux emp-single \
  --i-seqs emp-single-end-sequences.qza \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-column barcode-sequence \
  --o-per-sample-sequences demux.qza \
  --o-error-correction-details demux-details.qza
```

There are two main inputs for this process:
    
- `emp-single-end-sequences.qza`: Including fastq files of barcode and sequencing read.
- `sample-metadata.tsv`: Metatdata contains auxiliary information according to barcode. 

![input_files](static/input_files.png)

Base on the sequencing label on a read (in sequences fastq file), the original sample (in metadata) is retrieved by the barcode associating the same the sequencing label (in barcodes fastq file).  

## Demultiplex workflow

```{image} static/demux_workflow.png
:alt: demux_workflow
:height: 600px
:align: center
```

Sequencing reads are processed one by one from the input FASTQ file. For each read, the barcode sequence is first extracted from the beginning of the barcode read. If required, the barcode is reverse complemented to match the orientation of the reference barcodes.

Next, optional Golay error correction is applied to the barcode, allowing correction of sequencing errors (up to three mismatches) and improving robustness in sample identification. The corrected (or raw) barcode is then matched against a predefined barcode-to-sample mapping.

If a match is found, the corresponding read is assigned to that sample and written to its output FASTQ file. Reads that do not match any known barcode are discarded. Over the course of processing, reads are thus separated into multiple per-sample FASTQ files, each representing an individual sample from the original multiplexed dataset.