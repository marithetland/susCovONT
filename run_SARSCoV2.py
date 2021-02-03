#!/usr/bin/env python3
"""# marit.hetland@outlook.com
# github marithetland
# February 2021
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
    input_args.add_argument('-i', '--input_dir', type=pathlib.Path, required=True, help='Input directory, which should contain "fast5_pass" directory fast5-files.')
    #input_args.add_argument('-f', '--input_fast5', type=pathlib.Path, required=True, help='Input directory, which will be recursively searched for fast5-files.')
    #input_args.add_argument('-q', '--input_fastq', type=pathlib.Path, required=False, help='Input directory, which will be searched for fastq-files.')
    input_args.add_argument('-b', '--basecalling_model', type=str, required=True, choices=["r9.4_fast","r9.4_hac","r10_fast","r10_hac"], help='Indicate which basecalling mode to use. In most cases you probably want to use a HAC option.')
    input_args.add_argument('-k', '--barcode_kit', type=str, required=True, choices=["none","native_1-12","native_13-24","native_1-24"], help='Indicate which barcode-kits were used, if any.')
    input_args.add_argument('-s', '--sample_names', type=str, required=True, help='Provide a comma-separated list showing which barcode corresponds to which sample (for final report)')


    #Output - currently writes to same dir as input
    #output_args.add_argument('-o', '--outdir', type=pathlib.Path, required=False, default='.', help='Output directory for all output files. Only specify if different to input directory.')

    #Advanced options
    optional_args.add_argument('--resume_basecalling', action='store_true', required=False, help='Use this flag if your first run was interrupted and you want to resume. Default: off.')
    optional_args.add_argument('--cpu', action='store_true', required=False, help='If GPU is busy, use CPU with this flag. Will use 4 threads and 6 callers. Default: GPU.')
    #advanced_args.add_argument('-d', '--input_demultiplexed', type=pathlib.Path, required=True, help='Input directory of demultiplexed FASTQ files. Note: It is important that this demultiplexing was performed with the command described in the Wiki.')

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

#Only check if guppy_basecalling is being run
#def check_guppy_version():
#    run_command(['echo -n "guppy_basecaller \t" >> pipeline_versions.txt ; guppy_basecaller --version >> pipeline_versions.txt'], shell=True) ##Check for empty results, skip
#    pass

def check_nextflow_version():
    nextflow_version = run_command(['echo -n "nextflow \t" >> pipeline_versions.txt ; nextflow -version | grep version >> pipeline_versions.txt'], shell=True)     ##Check that version is correct
    current_version = "      version 20.10.0 build 5430"
    try:
        if nextflow_version == current_version:
            print("TADA")
    except:
        print("NOOO")
        sys.exit('Error: Nextflow version 20.10.0 is required')

def check_artic_version():
    run_command(['echo -n "artic \t" >> pipeline_versions.txt ;~/Programs/conda_for_covid/artic-2c6f8ebeb615d37ee3372e543ec21891/bin/artic --version >> pipeline_versions.txt'], shell=True) ##Check for empty results, skip
    pass

#def check_ncov2019_nextflow_exists():
#    if not os.path.exists(os.path.join(os.getcwd(),fast5_pass)):
#    pass


def fileCount(path, extension):
    count = 0
    for root, dirs, files in os.walk(path):
        count += sum(f.endswith(extension) for f in files)
    return count

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

    #Check that the input dir exists
    if not args.input_dir.is_dir():
        sys.exit('Error: {} is not a directory'.format(args.in_dir))
    #Check that the fast5_pass folder is present and count number of FAST5 files
    full_path = os.path.abspath(args.input_dir)
    fast5_pass=(full_path+"/fast5_pass/")
    fast5_pass_alt=(full_path+"/001_fast5_pass/")
    if not os.path.exists(os.path.join(os.getcwd(),fast5_pass)):
        if os.path.exists(os.path.join(os.getcwd(),fast5_pass_alt)):
            len_fast5=fileCount(fast5_pass_alt, '.fast5')
            if len_fast5 != 0:
                print("Found directory 001_fast5_pass with " + str(len_fast5) + " .fast5 files to analyse")
            else:
                sys.exit('Error: Found no .fast5 files')
        else:
            sys.exit('Error: {} is not a directory'.format(fast5_pass))
    if os.path.exists(os.path.join(os.getcwd(),fast5_pass)):
        len_fast5=fileCount(fast5_pass, '.fast5')
        if len_fast5 != 0:
            print("Found directory fast5_pass with " + str(len_fast5) + " .fast5 files to analyse")
        else:
            sys.exit('Error: Found no .fast5 files')

    #Check if the fastq_pass folder is present
    fastq_pass=(full_path+"/fastq_pass/")
    fastq_pass_alt=(full_path+"/002_basecalled/")
    if os.path.exists(os.path.join(os.getcwd(),fastq_pass)):
        print("Found fastq_pass folder. Checking for demultiplexed files (split into barcodes)")
    if not os.path.exists(os.path.join(os.getcwd(),fastq_pass)):
        if os.path.exists(os.path.join(os.getcwd(),fastq_pass_alt)):
            print("Found fastq_pass folder in directory 002_basecalled")
        else:
            print("No fastq_pass folder found. Will perform guppy basecalling and demultiplexing")

    #TODO: Check if Nextflow has already been run (007_nextflow + "nextflowSucess.txt")

    #TODO:Check barcode - sample doc

    #TODO: Add check for outdir if different than input dir


def get_guppy_basecalling_command(input_dir, save_dir, basecalling_model, resume, cpu, chunks):
    #guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i ../fast5_pass -s 002_basecalled -x auto -r 
    #Example command: guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i /path/to/reads -s run_name -x auto -r

    basecalling_command = ['guppy_basecaller ',
                     '--input_path ', input_dir, 
                     '--recursive ',
                     '--save_path ', save_dir]
    basecalling_command += BASECALLING[basecalling_model]
    if cpu:
        basecalling_command += ['--num_callers 6'] #change to auto / add option to use CPU or GPU with options
    else:
        basecalling_command += ['-x ', 'auto'] #change to auto / add option to use CPU or GPU with options

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

def get_nextflow_command(demultiplexed_fastq, fast5_pass, sequencing_summary,nf_outdir,run_name):
    #nextflow run ~/Programs/ncov2019-artic-nf/ -profile conda --prefix 20210202_1359_X5_FAO88697_5cf6e6f0  --cache /home/marit/Programs/conda_for_covid/work/conda --basecalled_fastq ./fastq_pass --fast5_pass  ./fast5_pass/ --sequencing_summary sequencing_summary_FAO88697_f6a6d889.txt --outdir 007_nextflow
    #nextflow run ~/Programs/ncov2019-artic-nf/ -profile conda --cache /home/marit/Programs/conda_for_covid/work/conda --prefix 210124_FAO88582_CoV_NB1-24   --basecalled_fastq ./003_demultiplexed --fast5_pass ./001_raw_fast5s/fast5_pass/new_copy     # --sequencing_summary ./002_basecalled/sequencing_summary.txt --outdir test_nextflow
    #TODO: add option to specify run folder, cache and medaka options
    nextflow_command = ['nextflow run ~/Programs/ncov2019-artic-nf -profile conda --nanopolish --cache /home/marit/Programs/conda_for_covid/work/conda ',
                     '--prefix', run_name, 
                     '--basecalled_fastq', demultiplexed_fastq, 
                     '--fast5_pass', fast5_pass,
                     '--sequencing_summary  ', sequencing_summary,
                     '--outdir ', nf_outdir]
    
    print(nextflow_command)
    return nextflow_command

def get_pangolin_command(consensus_dir):
    #conda update pangolin?
    #conda activate pangolin; pangolin file --outfile outfile ; pangolin deactivate
    #nextflow run ~/Programs/ncov2019-artic-nf/ -profile conda --cache /home/marit/Programs/conda_for_covid/work/conda --prefix 210124_FAO88582_CoV_NB1-24   --basecalled_fastq ./003_demultiplexed --fast5_pass ./001_raw_fast5s/fast5_pass/new_copy     # --sequencing_summary ./002_basecalled/sequencing_summary.txt --outdir test_nextflow
    pangolin_command = ['nextflow run ~/Programs/ncov2019-artic-nf -profile conda --nanopolish --cache /home/marit/Programs/conda_for_covid/work/conda ',
                     '--prefix', run_name, 
                     '--basecalled_fastq', demultiplexed_fastq, 
                     '--fast5_pass', fast5_pass,
                     '--sequencing_summary  ', sequencing_summary,
                     '--outdir ', nf_outdir]
    
    print(pangolin_command)
    return pangolin_command

def get_nextclade_command(consensus_dir):
    #docker pull image first
    #    #sudo docker run -it --rm -u 1001 --volume="${PWD}:/seq" neherlab/nextclade nextclade --input-fasta '/seq/sequences.fasta' --output-csv '/seq/results.csv'
    nextclade_command = ['nextflow run ~/Programs/ncov2019-artic-nf -profile conda --nanopolish --cache /home/marit/Programs/conda_for_covid/work/conda ',
                     '--prefix', run_name, 
                     '--basecalled_fastq', demultiplexed_fastq, 
                     '--fast5_pass', fast5_pass,
                     '--sequencing_summary  ', sequencing_summary,
                     '--outdir ', nf_outdir]
    
    print(nextclade_command)
    return nextclade_command

def yes_or_no(question):
    reply = str(input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        print("Starting pipeline")
        return 1
    elif reply[0] == 'n':
        print("Exiting pipeline")
        return sys.exit()
    else:
        return yes_or_no("Please Enter (y/n) ")

#main
def main():    
    args = parse_args()

    ## Set up log to stdout
    logfile= None
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m-%d-%Y %H:%M:%S')
    logging.info('Running covid-genomics v.1.0.0')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))


    ##Control input arguments
    full_path = os.path.abspath(args.input_dir) 
    run_name = os.path.basename(full_path) 
    
    outdir=full_path #TODO:Set optional (parental) outdir
    if outdir[-1] != '/':
        outdir = outdir + '/'

    #raw_fast5s=os.path.abspath(str(args.input_fast5))
    basecalled_fastq=(outdir+'002_basecalled/') #depends
    demultiplexed_fastq=(outdir+'003_demultiplexed/')
    nf_outdir=(outdir+'004_artic_minion/')
    consensus_dir=(outdir+'005_consensus_fasta/')
    sequencing_summary=(outdir+'002_basecalled/'+'sequencing_summary.txt') #depends on how it was done
    barcode_list=(args.sample_names)

    #Check that fast5_pass directory is present in input_dir:

    #TODO: Check if basecalling or basecalling+demultiplexing has been doneor not:



    #Check if it was done on GridION/miniT (i.e. fastq_pass/)

    #Check if it was done on Linux (i.e. 002_basecalled)

   #Specify barcode kits and basecaller mode
    barcode_kit=args.barcode_kit
    basecaller_mode=args.basecalling_model

    logging.info("##########CHECKPOINT##########")
    #Check necessary arguments, set variables and set up output directory
    check_python_version()
    check_guppy_version()
    check_nextflow_version()
    check_artic_version()
    #TODO: # check_ncov2019_nextflow_exists() #check only if not default

    ##Check with user that they are happy with the input for the pipeline
    print("Please check that the specified input below is correct for your run: ")
    check_arguments(args)
    print("The input folder is: " + full_path)
    print("The name of your run is: " + run_name)
    print("Your barcodes are specified in: " + barcode_list)
    print("Specified barcode kit is: " + barcode_kit)
    print("Using basecaller mode: " + basecaller_mode)
    while True:
        if(yes_or_no('Would you like to run this pipeline? y/n: ')):
            break
    ##

    ##Guppy basecalling
    #TODO: add check if basecalling has already been performed
    #if not args.basecalled:
        #Run basecalling
      #  basecalling_command=(get_guppy_basecalling_command(raw_fast5s, basecalled_fastq, basecaller_mode, args.resume_basecalling, args.cpu, args.chunks_per_runner))
     #   run_command([listToString(basecalling_command)], shell=True)
    #if not args.demultiplexed:
        ##Guppy demultiplexing
        #demultiplexing_command=(get_guppy_barcoder_command(basecalled_fastq, demultiplexed_fastq, barcode_kit))
        #run_command([listToString(demultiplexing_command)], shell=True)

    ##Run artic guppyplex and artic minion via PHW's nextflow pipeline
    #nextflow_command=(get_nextflow_command(demultiplexed_fastq, raw_fast5s, sequencing_summary,nf_outdir,run_name))
    #run_command([listToString(nextflow_command)], shell=True)

    ##Pangolin lineage assignment
    #subprocess.call(["activate", value])
    #subprocess.run('source activate environment-name && "enter command here" && source deactivate', shell=True)
    #consensus_dir
    #pangolin_command=(get_pangolin_command(consensus_dir))
    #run_command([listToString(pangolin_command)], shell=True)

    ##Nextclade lineage assignment and substitutions
    #As input: consensus_dir
    #nextclade_command=(get_nextclade_command(consensus_dir))
    #run_command([listToString(nextclade_command)], shell=True)
    #Make sure to pull image to download latest version
    #sudo docker run -it --rm -u 1001 --volume="${PWD}:/seq" neherlab/nextclade nextclade --input-fasta '/seq/sequences.fasta' --output-csv '/seq/results.csv'
    
    #QC: Generate report
    #with open(args.filename) as file:


if __name__ == '__main__':
    main()
#EOF
