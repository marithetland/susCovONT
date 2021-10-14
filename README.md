# susCovONT

**Pipeline to generate consensus.fasta files and identify pangolin lineage and nextstrain clade of Sars-CoV-2 genomes from ONT sequencing.**
* [Updates](#Updates)

This pipeline takes as input a folder with name `<run_name>` which contains the folders `fast5_pass` and `fastq_pass` and `sequencing_summary*.txt` from Sars-CoV-2 ONT sequencing together with a CSV-file which links barcode and sample name, and it [outputs consensus.fasta files](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) along with `<run_name>_report.csv` which includes [pangolin lineage](https://cov-lineages.org/pangolin.html), [nextstrain clade, mutations](https://clades.nextstrain.org/) and [QC](https://github.com/marithetland/susCovONT/wiki/3.-QC-and-parameters). 

## Installation
Install all necessary tools and environments with `susCovONT/scripts/install.sh`. You need to set the path of INSTALL_DIR, which is where the repositories will be installed (including this one), and [conda and docker has to be installed](https://github.com/marithetland/susCovONT/wiki/2.-Installation#installing-docker-and-conda):

```
INSTALL_DIR=/home/susamr/Programs/  #Change to your install dir
cd $INSTALL_DIR
git clone https://github.com/marithetland/susCovONT 
bash ./susCovONT/scripts/install.sh $INSTALL_DIR
```

## Basic usage

```
python susCovONT.py --input_dir /path/to/<run_name> --sample_names sample_names.csv
```
Where:
* `--input_dir`: Input directory `<run_name>` must contain `fast5_pass` and `fastq_pass` folders and `sequencing_summary*.txt`, with the `<run_name>` corresponding to your run (e.g. 20210213_1359_X5_FAO88697_5cf6e6f0)
* `--sample_names`: A CSV-file which connects barcodes with sample names, following the format:
```
barcode,sample_name
barcode01,NEGCONTROL
barcode02,E1234567_P1
NB03,V2345678_P1
```
The barcode column can take values following the format barcode[0-9][0-9] or NB[0-9][0-9] (as in the example above), and the sample_name column can be anything you'd like.


Note: Basecalling and demultiplexing may also be performed if not already done on GridION/MinIT.

## Please see the wiki for more information:
* [What does the pipeline do?](https://github.com/marithetland/covid-genomics/wiki/What-does-it-do%3F)
* [How to run](https://github.com/marithetland/covid-genomics/wiki/1.-How-to-run)
* [Installation](https://github.com/marithetland/covid-genomics/wiki/2.-Installation)
* [QC parameters and how to change these](https://github.com/marithetland/susCovONT/wiki/3.-QC-and-parameters)
* [What does the output look like](https://github.com/marithetland/covid-genomics/wiki/4.-Output)
* [How to run the commands manually](https://github.com/marithetland/covid-genomics/wiki/6.-Manual-run)

### And as an extra bonus:
* [Setting up MinKNOW and RAMPART](https://github.com/marithetland/covid-genomics/wiki/5.-MinKNOW-and-RAMPART)
* [Further exploration using the Nexctlade Web Application](https://github.com/marithetland/covid-genomics/wiki/Using-Nextclade-web-application)

## Note and thanks
Please note that this script was created for use at Stavanger University Hospital, you may need to change it (specifically the `scripts/config.cfg` file) for it to work in your environment.

This pipeline uses tools from the [Artic network's nCoV-2019 novel coronavirus bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html). See also the [QC and parameters](https://github.com/marithetland/susCovONT/wiki/3.-QC-and-parameters) page and [further links here](https://github.com/marithetland/covid-genomics/wiki/What-does-it-do%3F).

Many thanks to the artic, pangolin and nextclade developers for creating the protocols and pipelines!

## Repo name
This pipeline was created for the analysis of Sars-CoV-2 data from Oxford Nanopore Technologies (ONT) sequencing at Stavanger University Hospital (SUH/SUS). Hence the name, susCovONT: SUS + Covid-19 + ONT.

## Updates
- 2021-10-14: Added option for V3 or V4 primer schemes, default is now V4. Updated nextclade. Fixed output report.
- 2021-04-01: Combined the guppy basecalling and demultiplexing commands to match command (and output file structure) used by GridION
- 2021-03-31: Updated QC thresholds (see [QC and parameters](https://github.com/marithetland/susCovONT/wiki/3.-QC-and-parameters)) and added option for setting the `--normalise` value (default is `--normalise 500`) and also to automatically re-analyse samples that have 90-97% coverage without normalisng (`--renormalise=on`)
- 2021-03-12: Updated the threshold for QC PASS from >90% of bases confidently called (with 20X reads) to >97%, based on GISAID and [FHI's](https://github.com/folkehelseinstituttet/fhi-ncov-seq-pipelines) recommendations.
- ToDo: Merge guppy basecalling and demultiplexing commands to one to use same structure as ONT devices
- ToDo: Considering increasing the normalise value and/or using medaka in workflow instead of nanopolish
