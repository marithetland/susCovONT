# covid-genomics
Wrapper to basecall and demultiplex Sars-CoV-2 samples from ONT sequencing runs, then running the artic pipeline and finally assigning lineage with pangloin. 

This wrapper script is intended for the use at the AMR lab at Stavanger University Hospital. It performs basecalling and demultiplexing of fast5 files from a nanopore sequencing run. It takes as input the fast5-folder from your sequencing run and will output a consensus.fasta sequence which you can use for further analyses (e.g. phylogeny, pangolin, etc.).

At the bottom of this README file you will find the commands for each step of the script if you want to run the analyses outside the pipeline.

## Table of Contents
[Details](#Details)
[Computer requirements](#Requirements)  
[Input files](#Input)  
[Example Usage](#Basic-usage)  
[Usage](#Usage)  
[Output](#Output)

## Details
This script takes an ONT sequencing run, with or without barcodes, and performs basecalling and demultiplexing using the parameters described in https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html. It then runs the ARTIC pipeline via Nextflow (https://github.com/connor-lab/ncov2019-artic-nf) to produce consensus.fasta files for your Sars-CoV-2 genomes to be used in downstream analysis, such as pangolin lineage assignment or for upload to the GISAID database.

## QC
Enter details about QC.
20 reads?

## Requirements
These need to be installed and in path for the entire pipeline to work. Other versions of these tools will possibly work too, but these are the ones I have tested.

* guppy_basecaller (Available for customers at https://community.nanoporetech.com/)
* ncov2019-artic-nf (https://github.com/connor-lab/ncov2019-artic-nf)
Note: ncov2019-artic-nf currently comes with artic v1.1.3 as default. Change "artic=1.1.3" to "artic=1.2.1" in the file ncov2019-artic-nf/environments/nanopore/environment.yml to update to artic release v1.2.1.
For some reason the nextflow pipeline will not run unless the conda environment is created manually first. Do this by running:

```
conda env create --prefix /home/susamr/Programs/ncov2019-artic-nf/work7Conda/ --file ~/Programs/ncov2019-artic-nf/environments/nanopore/environment.yml 

```


## Input
* The script takes as input the folder where your fast5 files are placed. This can be in a folder called fast5_pass, or if using both fast5_pass and fast5_fail, use the parent directory as input, and it will recursively search for fast5 files.
* The script also needs to know where the ncov2019-artic-nf/ directory is placed. By default, the script assumes it is placed at ~/Programs/ncov2019-artic-nf/
* You must also specify barcodes: 1-12, 12-24 or both

## Usage:


## Output
The following directories are created:
* 002_basecalled - contains the basecalled FAST5 files (in FASTQ format)
* 003_demultiplexed - contains the demultiplexed FASTQ files (separated to barcodes)
* 004_artic-nextflow - contains the output of the artic + nextflow pipeline. See https://github.com/connor-lab/ncov2019-artic-nf and https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html for more information.



## Pangolin

Install pangolin as per the github instructions found here: https://github.com/cov-lineages/pangolin

```
cd ~/Programs
git clone https://github.com/cov-lineages/pangolin.git
cd pangolin
conda env create -f environment.yml
conda activate pangolin
python setup.py install
```


## Running it all manually
If for some reason the pipeline is not working or you want to run any of the steps manually, here are the steps to follow:

The code and text from steps 1-3 are from https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html

Assuming your fast5 reads are in a directory 001_raw_fast5:

1. Basecalling
```
guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i 001_raw_fast5 -s 002_basecalled -x auto -r
```

2. Demultiplexing

```
guppy_barcoder --require_barcodes_both_ends -i 002_basecalled -s 003_demultiplexed --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
```

3. Read filtering

```
artic guppyplex --min-length 400 --max-length 700 --directory 003_demultiplexed/barcode03 --prefix 004_guppyplex
```

This will perform a quality check. If you are only using “pass” reads you can speed up the process with:

```
conda activate artic
artic --version
artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory 003_demultiplexed/barcode03 --prefix 004_guppyplex/guppyplex
conda deactivate
```

4. Artic minion pipeline
```
conda activate artic
artic minion --normalise 200 --threads 4 --scheme-directory ~/Programs/artic-ncov2019/primer_schemes --read-file barcode03.fastq --fast5-directory 001_raw_fast5 --sequencing-summary sequencing_summary.txt nCoV-2019/V3 samplename
conda deactivate
```

5. Nextflow pipeline

The read filtering and the artic minion pipeline (steps 3 and 4) have to be run separately for each barcode. Matt Bull at PHW has created a nextflow pipeline which can take a whole sequencing run as input and runs these steps (3 and 4) on the entire run, leaving you with a consensus.fasta file for each barcode.

#Nextflow
```
nextflow run ~/Programs/ncov2019-artic-nf -profile conda --nanopolish --prefix nextflow --cache /media/susamr/lisa/SARS-CoV-2/test_nextflowArtic/work/conda/artic-2c6f8ebeb615d37ee3372e543ec21891/ --basecalled_fastq /media/susamr/lisa/SARS-CoV-2/sars-cov-2-analysis/nextflow_pipeline_check/002_basecalled_fast5_pass/ --fast5_pass /media/susamr/lisa/Nanopore_Runs/Cov-sekv/fast5/fast5_pass/ --sequencing_summary /media/susamr/lisa/SARS-CoV-2/sars-cov-2-analysis/nextflow_pipeline_check/002_basecalled_fast5_pass/sequencing_summary.txt
```

#TODO: Check the QC files


#Pangolin lineage assignment

When you are happy with your consensus.fasta file, you can run it through pangloin for lineage-assignment:


```
conda activate pangolin
pangolin ./005_/consensus_file.fasta --outfile ./006_pangolin_lineage/lineage_report.csv
conda deactivate pangolin
```

#END


## Updating tools
## TO DO:
Ideally, basecalling should be performed as sequencing is performed. Add option to skip basecalling if it has been done (or add option to add and set default to have been performed)
- Check for new releases of guppy, artic, the nextflow pipeline and pangolin every monday

To update pangolin:

```
conda activate pangolin
pangolin --update
```
##TODO: Should we use fast5_pass only or also fast5_fail?
##MAKE SURE VERSIONS ARE CHECKED AND LOGGED EVERY TIME IT IS RUN



#artic filtering in loop
for f in $(ls -d 003_demultiplexed/barcode[0-9][0-9]) ; do echo $f ; artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ${f} --prefix 004_guppyplex/guppyplex ; done


#artic minion pipeline in loop
artic minion --normalise 200 --threads 4 --scheme-directory ~/Programs/artic-ncov2019/primer_schemes --read-file ./004_guppyplex/test_guppyplex_barcode16.fastq --fast5-directory /media/susamr/maggie/ONT_covid/fast5_pass --sequencing-summary /media/susamr/maggie/ONT_covid/fastq_pass_old_v/MT-110273_20210124_182510_FAO88582_minion_sequencing_run_CoV_NB1-24_sequencing_summary.txt nCoV-2019/V3 ./005_consensus/barcode16

#Look into rampart



#In artic pipeline, add scheme and cache?
--schemeRepoURL /path/to/own/clone/of/github.com/artic-network/artic-ncov2019




##Input
Base 

#Filnavngivning
En bokstav (kode)


#ALSO adda a section about how to install all the tools

#ARTIC v1.2.1
#Nextflow pipeline, where:
Specify artic version
Install and then specify conda environment
Primer scheme too
