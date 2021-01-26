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

#import modules
import logging, time
import glob
import datetime
import itertools
import os, sys, re
import shutil
import csv
import pandas as pd
from argparse import ArgumentParser
import subprocess
from Bio import SeqIO
import pathlib
import collections
import dateutil.parser
import h5py
import random
import statistics
import tempfile
from subprocess import call
from subprocess import check_output, CalledProcessError, STDOUT

#Options for basecalling and barcoding ##Borrowed from Ryan Wick's  https://github.com/rrwick/MinION-desktop/blob/master/basecall.py
BASECALLING = collections.OrderedDict([
    ('r9.4_fast', ['--config dna_r9.4.1_450bps_fast.cfg ']), 
    ('r9.4_hac',  ['--config dna_r9.4.1_450bps_hac.cfg ']), #make default
    ('r10_fast',  ['--config dna_r10_450bps_fast.cfg ']),
    ('r10_hac',   ['--config dna_r10_450bps_hac.cfg ']),
])

BARCODING = collections.OrderedDict([
    ('native_1-12',  ['--arrangements_files "barcode_arrs_nb12.cfg "']),
    ('native_13-24', ['--arrangements_files "barcode_arrs_nb24.cfg" ']),
    ('native_1-24',  ['--arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" ']),
    #('rapid_1-12',   ['--barcode_kits "SQK-RBK004" --trim_barcodes ']),
    ('none',         [])
])

#Definitions
def parse_args():
    #Version
    parser = ArgumentParser(description='covid-genomics')
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + 'v.1.0.0')

    #Argsgroups
    input_args = parser.add_argument_group('Input options (required)')
    output_args = parser.add_argument_group('Output options (required)')
    optional_args = parser.add_argument_group('Optional flags')
    advanced_args = parser.add_argument_group('Advanced options')

    #Input
    input_args.add_argument('-i', '--input_fast5', type=pathlib.Path, required=True, help='Input directory, which will be recursively searched for fast5-files.')
    input_args.add_argument('-b', '--basecalling_model', type=str, required=True, choices=["r9.4_fast","r9.4_hac","r10_fast","r10_hac"], help='Indicate which basecalling mode to use. In most cases you probably want to use a HAC option.')
    input_args.add_argument('-k', '--barcode_kit', type=str, required=True, choices=["none","native_1-12","native_13-24","native_1-24"], help='Indicate which barcode-kits were used, if any.')

    #Output - currently writes to same dir as input
    output_args.add_argument('-o', '--outdir', type=pathlib.Path, required=False, default='.', help='Output directory for all output files. Only specify if different to input directory.')

    #Advanced options
    advanced_args.add_argument('--chunks_per_runner', type=str, required=False, help='Advanced option. Change chunks per runner. Default = 300')
    optional_args.add_argument('--resume', action='store_true', required=False, help='Use this flag if your first run was interrupted and you want to resume. Default: off.')
    optional_args.add_argument('--cpu', action='store_true', required=False, help='If GPU is busy, use CPU with this flag. Will use 4 threads and 6 callers. Default: GPU.')
    

    return parser.parse_args()

def listToString(s):  
    # initialize an empty string 
    str1 = " " 
    # return string   
    return (str1.join(s)) 


def run_command(command, **kwargs): #yekwahs
    command_str = ''.join(command)
    #logging.info('Running shell command: {}'.format(command_str))
    try:
        exit_status = call(command_str, **kwargs)
    except OSError as e:
        message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
        raise CommandError({"Error:": message})
    if exit_status != 0:
        message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
        raise CommandError({"Error:": message})


def check_python_version():
    try:
        assert sys.version_info >= (3, 5)
    except AssertionError:
        sys.exit('Error: Python 3.5 or greater is required')

def check_guppy_version():
    run_command(['guppy_basecaller --version > guppy_basecaller_version.txt'], shell=True) ##Check for empty results, skip
    pass

def check_nextflow_version():
    nextflow_version = run_command(['nextflow -version | grep version > nextflow_version.txt '], shell=True)     ##Check that version is correct
    current_version = "      version 20.10.0 build 5430"
    try:
        if nextflow_version == current_version:
            print("TADA")
    except:
        print("NOOO")
        sys.exit('Error: Nextflow version 20.10.0 is required')

def check_artic_version():
    run_command(['guppy_basecaller --version > guppy_basecaller_version.txt'], shell=True) ##Check for empty results, skip
    pass

def check_arguments(args):
    '''Check that the arguments are satisfied.'''
    barcode_choices = list(BARCODING.keys())
    args.barcode_kit = args.barcode_kit.lower()
    if args.barcode_kit not in barcode_choices:
        sys.exit('Error: valid --barcodes choices are: {}'.format(join_with_or(barcode_choices)))

    model_choices = list(BASECALLING.keys())
    args.basecalling_model = args.basecalling_model.lower()
    if args.basecalling_model not in model_choices:
        sys.exit('Error: valid --model choices are: {}'.format(join_with_or(model_choices)))

    if not args.input_dir.is_dir():
        sys.exit('Error: {} is not a directory'.format(args.in_dir))

    if args.outdir.is_file():
        sys.exit('Error: {} is a file (must be a directory)'.format(args.outdir))
    elif not args.outdir.is_file():
        directory=args.outdir
        if not os.path.exists(directory):
             os.makedirs(directory)
             print('Created output directory: '+ str(args.outdir))
        elif os.path.exists(args.outdir) and os.listdir(directory): 
            print("Warning: Specified output directory is not empty.")
        elif os.path.exists(args.outdir) and not os.listdir(directory):
            print("Specified output directory is empty.")
    
    return directory



def get_guppy_basecalling_command(input_dir, save_dir, basecalling_model, resume, cpu, chunks):
    #Example command: guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i /path/to/reads -s run_name -x auto -r

    basecalling_command = ['guppy_basecaller ',
                     '--input_path ', input_dir, 
                     '--recursive ',
                     '--save_path ', save_dir]
    basecalling_command += BASECALLING[basecalling_model]
    if cpu:
        basecalling_command += ['--num_callers 6'] #change to auto / add option to use CPU or GPU with options
    else:
        basecalling_command += ['--device ', 'cuda:all:100%'] #change to auto / add option to use CPU or GPU with options

    if chunks:
        basecalling_command += ['--chunks_per_runner ', chunks] #change to auto / add option to use CPU or GPU with options
    else:
        basecalling_command += ['--chunks_per_runner ', '300'] #change to auto / add option to use CPU or GPU

    if resume:
         basecalling_command += ['--resume']  #add option to resume
    
    print(basecalling_command)
    return basecalling_command

def get_guppy_barcoder_command(input_dir, save_dir, barcode_kit):
    #Example command: guppy_barcoder --require_barcodes_both_ends -i run_name -s output_directory --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"

    barcoding_command = ['guppy_barcoder ',
                     '--require_barcodes_both_ends '
                     '--input_path ', input_dir, 
                     '--save_path ', save_dir]
    barcoding_command += BARCODING[barcode_kit]
    #if resume:
    #     guppy_command += ['--resume']  #add option to resume
    
    print(barcoding_command)
    return barcoding_command

def get_nextflow_command(demultiplexed_fastq, fast5_pass, sequencing_summary):
    #nextflow run ~/Programs/ncov2019-artic-nf -profile conda --nanopolish --prefix nextflow --basecalled_fastq /media/susamr/lisa/SARS-CoV-2/sars-cov-2-analysis/nextflow_pipeline_check/002_basecalled_fast5_pass/ --fast5_pass /media/susamr/lisa/Nanopore_Runs/Cov-sekv/fast5/fast5_pass/ --sequencing_summary /media/susamr/lisa/SARS-CoV-2/sars-cov-2-analysis/nextflow_pipeline_check/002_basecalled_fast5_pass/sequencing_summary.txt

    nextflow_command = ['nextflow run ~/Programs/ncov2019-artic-nf -profile conda --nanopolish --prefix 004_artic-nf --cache  ',
                     '--basecalled_fastq', demultiplexed_fastq, 
                     '--fast5_pass', fast5_pass,
                     '--sequencing_summary  ', sequencing_summary]
    
    print(nextflow_command)
    return nextflow_command

#main
def main():    
    # Set up log to stdout
    logfile= None
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m-%d-%Y %H:%M:%S')
    logging.info('Running susamr-bascealling script v.1.0.0')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))

    #Check arguments, set variables and set up output directory

    logging.info("Hello! Just checking that all parameters and versions are in place before we start...")
    args = parse_args()
    check_python_version()
    check_guppy_version()
    check_nextflow_version()
    #check_arguments()
    #TODO: MUST also check artic version and ncov-nf folder is correct or exists:
    #check_ncov2019_nextflow_exists() #check only if not default
    #check_artic_version()
    #Must also set outdir as optional ...

    outdir=os.path.abspath(str(args.outdir))
    if outdir[-1] != '/':
        outdir = outdir + '/'
        print(outdir)

    #Specify output directories (this can be done better)
    raw_fast5s=os.path.abspath(str(args.input_fast5))
    basecalled_fastq=(outdir+'002_basecalled/')
    demultiplexed_fastq=(outdir+'003_demultiplexed/')
    #prefix=(outdir+'004_artic-nf/')

    barcode_kit=args.barcode_kit
    print("Specified barcode kit is: " + barcode_kit)
    basecaller_mode=args.basecalling_model
    print("Using basecaller mode: " + basecaller_mode)



    
    logging.info("We're good to go! Basecalling and demultiplexing with Guppy now")
    logging.info("PS. You can see which version of guppy was used in the file: guppy_basecaller_version.txt")

    ##Part 1: Run Guppy    

    ##Guppy basecalling and demultiplexing
    ##Is performed as described on: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html

    ##Command:
    #guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i /path/to/reads -s run_name -x auto -r
    basecalling_command=(get_guppy_basecalling_command(raw_fast5s, basecalled_fastq, basecaller_mode, args.resume, args.cpu, args.chunks_per_runner))
    #run_command([listToString(basecalling_command)], shell=True)
    ##Guppy demultiplexing
    #Command
    #guppy_barcoder --require_barcodes_both_ends -i run_name -s output_directory --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
    #get_guppy_barcoder_command()
    demultiplexing_command=(get_guppy_barcoder_command(basecalled_fastq, demultiplexed_fastq, barcode_kit))
    #run_command([listToString(demultiplexing_command)], shell=True)

    ##Nextflow pipeline
    ##Notes:
    #command
    #nextflow run ~/Programs/ncov2019-artic-nf -profile conda --nanopolish --prefix nextflow --cache /media/susamr/lisa/SARS-CoV-2/test_nextflowArtic/work/conda/artic-2c6f8ebeb615d37ee3372e543ec21891/ --basecalled_fastq /media/susamr/lisa/SARS-CoV-2/sars-cov-2-analysis/nextflow_pipeline_check/002_basecalled_fast5_pass/ --fast5_pass /media/susamr/lisa/Nanopore_Runs/Cov-sekv/fast5/fast5_pass/ --sequencing_summary /media/susamr/lisa/SARS-CoV-2/sars-cov-2-analysis/nextflow_pipeline_check/002_basecalled_fast5_pass/sequencing_summary.txt
    sequencing_summary=(outdir+'002_basecalled/'+'sequencing_summary.txt')

    nextflow_command=(get_nextflow_command(demultiplexed_fastq, raw_fast5s, sequencing_summary))
    run_command([listToString(nextflow_command)], shell=True)


    ##Pangolin lineage assignment
    #subprocess.call(["activate", value])
    #subprocess.run('source activate environment-name && "enter command here" && source deactivate', shell=True)




if __name__ == '__main__':
    main()
#EOF
