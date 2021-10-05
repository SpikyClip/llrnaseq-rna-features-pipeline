# llrnaseq-rna-features-pipeline

This readme explains how to use the Nextflow llrnaseq in conjunction
with the rna-features python package to generate transfer learning
expression features.

## Process
1. Install [`Nextflow`](https://nf-co.re/usage/installation)
   (`>=21.04.3`) and the
   [`llrnaseq`](https://github.com/SpikyClip/llrnaseq#quick-start)
   pipeline. 
2. Install [rna-features](https://github.com/SpikyClip/rna-features#installation)
   (`python>=3.9`).
3. Create a [`samplesheet.csv`](https://github.com/SpikyClip/llrnaseq/blob/master/docs/usage.md#samplesheet-input)
   that points to your `fastq.gz reads`.
4. Run `llrnaseq` with the appropriate
   [options](https://github.com/SpikyClip/llrnaseq#quick-start) for the
   sample species.
5. Perform the appropriate differential expression analysis on the gene
   `counts.txt` produced by `llrnaseq` using
   [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
   A sample R script along with inputs and expected output contrast
   files can be found in the `deseq2_example` folder.
6. Repeat steps 3 - 5 on each dataset, storing the `tpm.tsv` (produced
   by `llrnaseq`) and contrast files (produced by `DESeq2`) in folders by dataset.
7. [Generate an expression features matrix](https://github.com/SpikyClip/rna-features#usage) by passing the dataset
   directory paths to rna-features.
