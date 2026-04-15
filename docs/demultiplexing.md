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

Sequencing reads are processed one by one from the input FASTQ file. For each read, the associated barcode sequence is retrieved. Although, the step "Reverse complement barcode" was not performed in this command (`demux emp-single`), more information can be found in [](#reverse-complement).

Next, Golay error correction is applied by default to the barcode, allowing correction of sequencing errors in barcode sequence (up to three mismatches) and improving robustness in sample identification (read more in [](#golay-correct)). 

The corrected barcode is then matched against a predefined barcode-to-sample mapping. If a match is found, the corresponding read is assigned to that sample and written to its output FASTQ file. Reads that do not match any known barcode are discarded. 

Over the course of processing, reads are thus separated into multiple per-sample FASTQ files, each representing an individual sample from the original multiplexed dataset.

(reverse-complement)=
## Reverse complement barcode

(golay-correct)=
## Golay error correction

Golay error correction is an error-correcting coding method that represents DNA barcodes as structured 24-bit vectors, enabling detection and correction of sequencing errors based on the properties of the Golay (24,12,8) code [@morelos-zaragozaArtErrorCorrecting2006, pg.30-31].

It is used in sequencing workflows to improve sample assignment accuracy. Because sequencing errors are common, exact barcode matching can lead to substantial data loss. Golay coding allows correction of up to three bit errors, thereby recovering reads while minimizing incorrect assignments.

**Correction workflow:**
1. Convert barcode (DNA) to a 24-bit vector (2 bits per nucleotide)  
2. Compute the syndrome using a parity-check matrix \(H\)  
3. Use the syndrome to identify the most likely error pattern (lookup table)  
4. Correct the bit vector by applying the error pattern  
5. Convert the corrected bits back to a DNA sequence  
6. Match the corrected barcode to a known sample identifier  

Golay error correction assumes that true barcodes belong to the predefined set of valid Golay codewords, that sequencing errors are limited (typically no more than three bit errors), and that these errors occur randomly rather than in a systematic or biased manner.

Full python script for Golay error correction can be found in [GolayDecoder.py](https://github.com/Truongphi20/qiime2_blog/blob/demultiplex/algorithm_reference/GolayDecoder.py).