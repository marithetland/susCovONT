#!/usr/bin/env python3
"""# marit.hetland@outlook.com
# github marithetland
# January 2021
# This script will basecall and demultiplex
# fast5 files from an ONT sequencing run before
# running the ARTIC minion pipeline via PHW's
# Nextflow pipeline. 
# Currently tested with versions
# Guppy basecalling software: Version 4.2.2+effbafg
# Nextflow v.20.10.0
# Artic v1.2.1
# ncov nextflow pipeline last updated 23.12.2020 (https://github.com/connor-lab/ncov2019-artic-nf)
"""

#####################################
### Sars-CoV-2 genomics @SUS
### Marit Hetland, January 2021
###
### This script will first guppy basecall and demultiplex fast5 reads from a nanopore run
### It will then run the sequencing run through the Nextflow artic pipeline which currently uses
### the artic minion pipeline v.1.2.1 
### Nextflow version = 20.10.0
#####################################


##Guppy basecalling and demultiplexing
##Is performed as described on: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html

##Guppy basecalling - note we use high accuracy base calling

##Command:
#guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i /path/to/reads -s run_name -x auto -r

##Guppy demultiplexing
#Command
#guppy_barcoder --require_barcodes_both_ends -i run_name -s output_directory --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"

##Nextflow pipeline
##Notes:
#command
#nextflow run ~/Programs/ncov2019-artic-nf -profile conda --nanopolish --prefix nextflow --basecalled_fastq /media/susamr/lisa/SARS-CoV-2/sars-cov-2-analysis/nextflow_pipeline_check/002_basecalled_fast5_pass/ --fast5_pass /media/susamr/lisa/Nanopore_Runs/Cov-sekv/fast5/fast5_pass/ --sequencing_summary /media/susamr/lisa/SARS-CoV-2/sars-cov-2-analysis/nextflow_pipeline_check/002_basecalled_fast5_pass/sequencing_summary.txt





