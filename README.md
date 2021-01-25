# covid-genomics
Wrapper to basecall, demultiplex and run the artic pipeline for Sars-CoV-2 sequences

This wrapper script is intended for the use at the AMR lab at Stavanger University Hospital. It performs basecalling and demultiplexing of fast5 files from a nanopore sequencing run. It takes as input the fast5-folder from your sequencing run and will output a consensus.fasta sequence which you can use for further analyses (e.g. phylogeny, pangolin, etc.).


## Table of Contents
[Details](#Details)
[Computer requirements](#Requirements)  
[Input files](#Input-files)  
[Example Usage](#Basic-usage)  
[Usage](#Usage)  
[Output](#Output)

## Details
This script takes an ONT sequencing run, with or without barcodes,

## Requirements
These need to be installed and in path for the entire pipeline to work. Other versions of these tools will possibly work too, but these are the ones I have tested.

* guppy_basecaller (Available for customers at https://community.nanoporetech.com/)
* ncov2019-artic-nf (https://github.com/connor-lab/ncov2019-artic-nf)
Note: ncov2019-artic-nf currently comes with artic v1.1.3 as default. Change "artic=1.1.3" to "artic=1.2.1" in the file ncov2019-artic-nf/environments/nanopore/environment.yml to update to artic release v1.2.1.
For some reason the nextflow pipeline will not run unless the conda environment is created manually first. Do this by running:

```
conda env create --prefix /home/susamr/Programs/ncov2019-artic-nf/work7Conda/ --file ~/Programs/ncov2019-artic-nf/environments/nanopore/environment.yml 

```


## Input files
* The script takes as input the folder where your fast5 files are placed. This can be in a folder called fast5_pass, or if using both fast5_pass and fast5_fail, use the parent directory as input, and it will recursively search for fast5 files.
* The script also needs to know where the ncov2019-artic-nf/ directory is placed. By default, the script assumes it is placed at ~/Programs/ncov2019-artic-nf/

## Usage:


## Output
The following directories are created:
* 002_basecalled - contains the basecalled FAST5 files (in FASTQ format)
* 003_demultiplexed - contains the demultiplexed FASTQ files (separated to barcodes)
* 004_artic-nextflow - contains the output of the artic + nextflow pipeline. See https://github.com/connor-lab/ncov2019-artic-nf and https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html for more information.
