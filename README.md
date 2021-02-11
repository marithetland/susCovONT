# covid-genomics

This pipeline takes as input Sars-CoV-2 ONT sequencing data and outputs consensus.fasta files along with assigned lineage ([pangolin](https://cov-lineages.org/pangolin.html) and [nextclade](https://clades.nextstrain.org/)), mutations and QC. It uses scripts and pipelines from the [Artic network's nCoV-2019 novel coronavirus bioinformatics protocol] (https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

Takes as input a directory containing fast5_pass with/without fastq_pass directories + a file to link sample names and barcodes.

Please see the wiki for how to:
## Table of Contents
* [What does the pipeline do?](https://github.com/marithetland/covid-genomics/wiki/What-does-it-do%3F)
* [Install this program and its dependencies](https://github.com/marithetland/covid-genomics/wiki/Installation)
* [How to run the pipeline and required inputs](https://github.com/marithetland/covid-genomics/wiki/How-to-run)
* [The QC parameters used and how to change these](https://github.com/marithetland/covid-genomics/wiki/QC)
* [Example ouput report and explanation](https://github.com/marithetland/covid-genomics/wiki/Output)

Additionally, you might be interested in:
* Effectively copy files from GridION/MinIT to analysis computer
* [How to run each command manually](https://github.com/marithetland/covid-genomics/wiki/Manual-run)

Please note that this script is intended for the use at the AMR lab at Stavanger University Hospital, you may need to change it for it to work in your environment.

Many thanks to the Artic network team, pangolin and nextclade developers for developing the protocols and pipelines used for this!
