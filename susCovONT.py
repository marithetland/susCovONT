#!/usr/bin/env python3
""" marit.hetland@outlook.com
github marithetland, February 2021
This is a pipeline to generate consensus.fasta files 
and to identify pangolin lineage and nextstrain clade of 
Sars-CoV-2 genomes from ONT sequencing.
Can also perform guppy basecalling and demultiplexing.
It uses the ARTIC minion pipeline via PHW's Nextflow pipeline. 
Currently tested with versions
- Guppy basecalling software: Version 4.4.1
- Nextflow v.20.10.0
- Artic v1.2.1
- ncov nextflow pipeline last updated 23.12.2020 (https://github.com/connor-lab/ncov2019-artic-nf)
"""

#import modules
import logging, time
import glob
import datetime
import itertools
import os, sys, re
import configparser
import shutil
import csv
import pandas as pd
from functools import reduce
from argparse import ArgumentParser
import subprocess
#from Bio import SeqIO
import pathlib
import collections
#import dateutil.parser
#import random
#import statistics
#import tempfile
from subprocess import call, check_output, CalledProcessError, STDOUT
from pathlib import Path    
from shutil import copyfile

#Options for basecalling and barcoding
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
    ('native_1-96',  ['--arrangements_files "barcode_arrs_nb96.cfg" ']),
    #('rapid_1-12',   ['--barcode_kit "SQK-RBK004" --trim_barcodes ']),
    ('none',         [])
])

#Definitions
def parse_args():
    #Version
    parser = ArgumentParser(description='susCovONT')
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + 'v.1.0.0')

    #Argsgroups
    input_args = parser.add_argument_group('Input options (required)')
    advanced_args = parser.add_argument_group('Advanced options')
    optional_args = parser.add_argument_group('Necessary flags only if performing guppy basecalling/demultiplexing')
    basecalling_args = parser.add_argument_group('Options for basecalling command')
    #output_args = parser.add_argument_group('Output options (required)')

    #Input
    input_args.add_argument('-i', '--input_dir', type=pathlib.Path, required=True, help='Input directory, which should contain "fast5_pass" directory fast5-files.')
    input_args.add_argument('-s', '--sample_names', type=pathlib.Path, required=True, help='Provide a comma-separated list showing which barcode corresponds to which sample (for final report)')
    
    #Advanced options
    advanced_args.add_argument('--generate_report_only', action='store_true', required=False, help='Do not run any tools, just (re)generate output report from already completed run. Default: off.')
    advanced_args.add_argument('--offline', action='store_true', required=False, help='The script downloads the newest primer schemes, nextclade and pangolin each time it runs. Use this flag if you want to run offline with already installed versions.fault: off.')
    advanced_args.add_argument('--no_move_files', action='store_true', required=False, help='By default, the input fast5_pass and fastq_pass dirs will be moved to subdir 001_rawData. Use this flag if you do not want that')
    advanced_args.add_argument('--no_artic', action='store_true', required=False, help='Use this flag to run only pangolin and nextclade on an already completed artic nextflow (with same folder structure)')
    advanced_args.add_argument('--dry_run', action='store_true', required=False, help='Executes nothing. Prints the commands that would have been run in a non-dry run.')

    #Options to perform basecalling and/or demultiplexing before running artic guppyplex and minion
    optional_args.add_argument('-b', '--basecalling_model', type=str, required=False, choices=["r9.4_fast","r9.4_hac","r10_fast","r10_hac"], help='Use flag to perform basecalling before running the artic pipeline. Indicate which basecalling mode to use. In most cases you want to use a HAC option.')
    optional_args.add_argument('-k', '--barcode_kit', type=str, required=False, choices=["none","native_1-12","native_13-24","native_1-24","native_1-96"], help='Use flag to perform demultiplexing of basecalled data. Indicate which barcode-kits were used (or none).')
    basecalling_args.add_argument('--guppy_resume_basecalling', action='store_true', required=False, help='This flag can be used with --basecalling to resume an interrupted basecalling run. Default: off.')
    basecalling_args.add_argument('--guppy_use_cpu', action='store_true', required=False, help='This flag can be used with --basecalling to run on CPU instead of GPU. Will use 4 threads and 6 callers. Default: GPU -auto x.')

    #Output - currently writes to same dir as input
    #output_args.add_argument('-o', '--outdir', type=pathlib.Path, required=False, default='.', help='Output directory for all output files. Only specify if different to input directory.')
    return parser.parse_args()

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

def listToString(s):  
    # initialize an empty string 
    str1 = " " 
    # return string   
    return (str1.join(s)) 

def combineCommand(s):  
    # initialize an empty string 
    str2 = "" 
    # return string   
    return (str2.join(s)) 

def join_with_or(str_list):
    if isinstance(str_list, dict):
        str_list = list(str_list.keys())
    if len(str_list) == 0:
        return ''
    if len(str_list) == 1:
        return str_list[0]
    return ', '.join(str_list[:-1]) + ' or ' + str_list[-1]

def fileCount(path, extension):
    count = 0
    for root, dirs, files in os.walk(path):
        count += sum(f.endswith(extension) for f in files)
    return count

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

def check_python_version():
    try:
        assert sys.version_info >= (3, 5)
        run_command(['echo "Pipeline versions \n" > pipeline_versions.txt'], shell=True) #Overwrite any previous file
    except AssertionError:
        sys.exit('Error: Python 3.5 or greater is required')

def check_guppy_version():
   run_command(['echo -n "guppy_basecaller \t" >> pipeline_versions.txt ; guppy_basecaller --version >> pipeline_versions.txt'], shell=True) ##Check for empty results, skip
   pass

def check_nextflow_version():
    nextflow_version = run_command(['echo -n "nextflow \t" >> pipeline_versions.txt ; nextflow -version | grep version >> pipeline_versions.txt'], shell=True)     ##Check that version is correct
    current_version = "      version 20.10.0 build 5430"
    try:
        if nextflow_version == current_version:
            print("Nextflow version is acceptable")
    except:
        sys.exit('Error: Nextflow version 20.10.0 is required')
    return

def check_artic_version(conda_location):
    run_command(['echo -n "artic \t" >> pipeline_versions.txt ; ',conda_location,'/artic-2c6f8ebeb615d37ee3372e543ec21891/bin/artic --version >> pipeline_versions.txt'], shell=True) ##Check for empty results, skip
    pass

# def check_pangolin_version():
#     #TODO
#     pass

# def check_ncov2019_nextflow_exists(nf_dir_location):
#     #TODO
#     nf_dir_path = Path(nf_dir_location)
#     if not nf_dir_path_full.is_dir():
#         sys.exit('Error: {} is not a directory'.format(nf_dir_path))

def set_config_variables():
    base_dir = os.path.dirname(os.path.realpath(__file__))
    #Sort out variables from config file
    config=configparser.ConfigParser()
    scripts_path=(base_dir+'/scripts')
    config.read(os.path.join(scripts_path,'config.cfg'))
    config.sections()
    nf_dir_location=config['NF_LOCATION']['nf_location']
    conda_location=config['CONDA_LOCATION']['conda_location']
    number_CPUs=config['CPU']['cpus']
    schemeRepoURL=config['OFFLINE']['schemeRepoURL']

    #Set CPU for nextflow
    resources_file=(nf_dir_location+'conf/resources.config')
    set_cpu_to=str("        cpus = "+number_CPUs)
    set_cpu_command = ["sed -i '5s/.*/",set_cpu_to,
                       "/' ",resources_file]
    run_command([combineCommand(set_cpu_command)], shell=True)

    return nf_dir_location, conda_location, number_CPUs, schemeRepoURL


def check_what_to_run(args):
    if args.input_dir:
        full_path = os.path.abspath(args.input_dir)
        run_name = os.path.basename(full_path) 
        outdir=full_path #TODO:Set optional (parental) outdir
        if outdir[-1] != '/':
            outdir = outdir + '/'
        nf_outdir=(outdir+'002_articPipeline/')
        artic_outdir=(outdir+'002_articPipeline/qc_pass_climb_upload/')
        consensus_dir=(outdir+'003_consensusFasta/')
        pangolin_outdir=(outdir+'004_pangolin/')
        nextclade_outdir=(outdir+'005_nextclade/')
        nextclade_outfile=(outdir+'005_nextclade/'+run_name+'_sequences.fasta_nextclade.csv')
        pangolin_outfile=(outdir+'004_pangolin/lineage_report.csv')
        artic_qc=(outdir+'002_articPipeline/'+run_name+'.qc.csv')
        sample_name=str(args.sample_names)
        raw_data_path=(outdir+'001_rawData/')
        final_report_name=(run_name+'_report.csv')
    else:
        sys.exit("Error: No input folder specified.")

    return full_path, run_name, outdir, nf_outdir, artic_outdir, consensus_dir, pangolin_outdir, nextclade_outdir, nextclade_outfile, pangolin_outfile, artic_qc,sample_name,final_report_name, raw_data_path

def check_arguments(args):
    ''' Check which arguments are specified and apply '''
    ##First check for presence of input dir with fast5 and (demultiplexed) fastq files
    fast5_pass_path, fastq_pass_path, fastq_pass_dem_path, sequencing_summary_path, print_len_fast5 = check_input(args)

    ##Config
    nf_dir_location, conda_location, number_CPUs, schemeRepoURL = set_config_variables()

    ##Check that tools are installed, set variables
    check_python_version()
    check_nextflow_version()
    check_artic_version(conda_location) #get from config
    #TODO check_ncov2019_nextflow_exists(nf_dir_location) #get from config

    ##TODO: Set up optional output directory

    ##Check that basecalling/demultiplexing command have been properly specified
    basecalling_model, barcode_kit = check_guppy_params(args)

    ##TODO:Check that the sample_names.csv file exists and follows correct format
    
    ##Create output directories and variables
    full_path, run_name, outdir, nf_outdir, artic_outdir, consensus_dir, pangolin_outdir, nextclade_outdir, nextclade_outfile, pangolin_outfile, artic_qc, sample_name, final_report_name, raw_data_path = check_what_to_run(args)

    logging.info("##########CHECKPOINT##########")
    print("The input folder is: " + full_path)
    print("The name of your run is: " + run_name)
    print("Found fast5_pass directory with "+print_len_fast5+" fast5 files to analyse")
    pipeline_commmand = ['The susCovONT pipeline will: ']
    if args.basecalling_model:
        pipeline_commmand += ['\n- run basecalling with model', basecalling_model] 
    if args.barcode_kit:
        pipeline_commmand += ['\n- run demultiplexing with barcodes ', barcode_kit] 
    if not args.no_artic:
        pipeline_commmand += ['\n- run the artic pipeline, pangolin and nextclade '] 
    if args.no_artic:
        pipeline_commmand += ['\n- run pangolin and nextclade on already completed artic run '] #TODO add check for existing run
    if args.generate_report_only:
        pipeline_commmand += ['\n- only (re)create output report from already completed pipeline run'] 
    if args.offline:
        pipeline_commmand += ['\n- in offline mode'] 
    if args.dry_run:
        pipeline_commmand += ['\n- in dry mode (no execution of commands)'] 

    print(listToString(pipeline_commmand))

    while True:
        if(yes_or_no('Would you like to run this pipeline? y/n: ')):
            break

    return full_path, fast5_pass_path, fastq_pass_path, fastq_pass_dem_path, sequencing_summary_path, basecalling_model, barcode_kit, run_name, outdir, nf_outdir, artic_outdir, consensus_dir, pangolin_outdir, nextclade_outdir, nextclade_outfile, pangolin_outfile, artic_qc, sample_name, final_report_name, raw_data_path

def get_sample_names(args):
    if not args.sample_names.is_file():
        sys.exit('Error: {} is not a file. Please check your input.'.format(args.sample_names))
    if args.sample_names.is_file():
        sample_names=str(args.sample_names) #TODO: Must check that file exists and has correct header
        print("Your barcodes are specified in: " + os.path.abspath(sample_names))
    
    #TODO: Check that it is the correct format
    return sample_names

def check_input(args): #TODO: This can be done neater
    ''' Directory with FAST5 files must be present for thepipeline to run '''
    #FAST5 files must always be present. Look recursively for FAST5 files.
    #Check that the input dir exists
    if not args.input_dir.is_dir():
        sys.exit('Error: {} is not a directory'.format(args.in_dir))
    
    ##Check FAST5
    #Check that the fast5_pass folder is present and count number of FAST5 files
    full_path = os.path.abspath(args.input_dir)
    fast5_pass=(full_path+"/fast5_pass/") 
    fast5_pass_alt=(full_path+"/001_rawData/fast5_pass/")
    if not os.path.exists(os.path.join(os.getcwd(),fast5_pass)):
        if os.path.exists(os.path.join(os.getcwd(),fast5_pass_alt)):
            len_fast5=fileCount(fast5_pass_alt, '.fast5')
            if len_fast5 != 0:
                print_len_fast5=str(len_fast5)
                fast5_pass_path=fast5_pass_alt
            else:
                sys.exit('Error: Found no .fast5 files. Please check that your input directory is correct.')
        else:
            sys.exit('Error: {} is not a directory'.format(fast5_pass))
    if os.path.exists(os.path.join(os.getcwd(),fast5_pass)):
        len_fast5=fileCount(fast5_pass, '.fast5')
        if len_fast5 != 0:
            print_len_fast5=str(len_fast5)
            fast5_pass_path=fast5_pass
        else:
            sys.exit('Error: Found no .fast5 files')

    ''' Unless basecalling or demultiplexing is being performed, FASTQ files
    must be present for the pipeline to run '''
    #Check if the fastq_pass folder is present in either of:
    fastq_pass=(full_path+"/fastq_pass") #FASTQ demultiplexed on ONT (guppy)
    fastq_pass_undem=(full_path+"/001_rawData/fastq_pass_notDemultiplexed") #FASTQ PASS from guppy basecaller, not demultiplexed. Exists together with fastq_pass
    fastq_pass_dem=(full_path+"/001_rawData/fastq_pass") #FASTQ PASS from guppy basecaller, demultiplexed
    sequencing_summary=""

    ##IF not demultiplexing or basecalling:
    if not args.barcode_kit and not args.basecalling_model:
        #The fastq_pass folder must exist as we are running nextflow directly on this + fast5.
        if os.path.exists(os.path.join(os.getcwd(),fastq_pass_dem)):
            fastq_pass_dem_path=fastq_pass_dem
            if os.path.exists(os.path.join(os.getcwd(),fastq_pass_undem)):
                sequencing_summary=(fastq_pass_undem+"/sequencing_summary.txt") 
                fastq_pass_path=fastq_pass_undem
            elif not os.path.exists(os.path.join(os.getcwd(),fastq_pass_undem)):
                sequencing_summary=(full_path+"/sequencing_summary*.txt") 
                fastq_pass_path=fastq_pass_dem
            else:
                sys.exit('Error: Could not find sequencing_summary.txt file to match ' + fastq_pass_dem)
        elif os.path.exists(os.path.join(os.getcwd(),fastq_pass)):
            #If already basecalled AND demultiplexed on the GridION/MinIT:
            fastq_pass_path=fastq_pass
            fastq_pass_dem_path=fastq_pass
            sequencing_summary=(full_path+"/sequencing_summary*.txt")  #TODO: split run name to the last two "_"s to get this name properly
    ##IF demultiplexing
    if args.barcode_kit:
        #Check that there isn't already a folder called fastq_pass_demultiplexed
        if os.path.exists(os.path.join(os.getcwd(),fastq_pass_dem)):
            sys.exit('Error: The folder {} already exists. Please delete/move this folder if you want to perform demultiplexing.'.format(str(fastq_pass_dem)))
    ##IF basecalling and demultiplexing
    if args.barcode_kit and args.basecalling_model:
        if os.path.exists(os.path.join(os.getcwd(),fastq_pass)):
            sys.exit('Error: The folder {} already exists. Please delete/move this folder if you want to perform basecalling.'.format(str(fastq_pass)))
        elif os.path.exists(os.path.join(os.getcwd(),fastq_pass_dem)):
            sys.exit('Error: The folder {} already exists. Please delete/move this folder if you want to perform basecalling.'.format(str(fastq_pass_dem)))
        elif os.path.exists(os.path.join(os.getcwd(),fastq_pass_undem)):
            sys.exit('Error: The folder {} already exists. Please delete/move this folder if you want to perform basecalling.'.format(str(fastq_pass_undem)))
        else:
            fastq_pass_path=fastq_pass_undem
            fastq_pass_dem_path=fastq_pass_dem
            sequencing_summary=(fastq_pass_undem+"/sequencing_summary.txt") 
    ##IF demultiplexing but NOT basecalling
    if args.barcode_kit and not args.basecalling_model:
        #If basecalling not specified, then the basecalled fastq folder must already exist
        if os.path.exists(os.path.join(os.getcwd(),fastq_pass_undem)):
            fastq_pass_path=fastq_pass_undem
            sequencing_summary=(fastq_pass_undem+"/sequencing_summary.txt")
            fastq_pass_dem_path=fastq_pass_dem
        elif os.path.exists(os.path.join(os.getcwd(),fastq_pass)):
            fastq_pass_dem_path=fastq_pass
            sequencing_summary=(full_path+"/sequencing_summary*.txt")  #TODO: split run name to the last two "_"s to get this name properly
        else:
            sys.exit("Error: Did not find a fastq_pass folder in {}/fastq_pass or {}/001_rawData/fastq_pass. Have you basecalled your fast5s or did you forget to specify basecalling (--basecalling_model)?").format(str(input_dir))
    ##IF demultiplexing but NOT basecalling specified
    if args.basecalling_model and not args.barcode_kit:
        sys.exit("Error: Cannot perform basecalling and downstream artic analysis without demultiplexing. Please also specify --barcode_kit /-k in your command.")

    #Check that the specified sequencing summary actually exists
    if not args.basecalling_model or not args.barcode_kit:
        seq_sum_path=Path(sequencing_summary)
        if not glob.glob(sequencing_summary):
            sys.exit('Error: {} is not a file'.format(sequencing_summary))


    return fast5_pass_path, fastq_pass_path, fastq_pass_dem_path, sequencing_summary, print_len_fast5

def check_guppy_params(args): 
    model_choices = list(BASECALLING.keys())
    barcode_choices = list(BARCODING.keys())
    basecalling_model=""
    barcode_kit=""
    #Check if basecalling is needed and if correct parameters are given
    if args.basecalling_model:
        basecalling_model = args.basecalling_model.lower()
        check_guppy_version()
        if args.basecalling_model not in model_choices:
            sys.exit('Error: valid --model choices are: {}'.format(join_with_or(model_choices)))

    #Check if demultiplexing is needed and if correct parameters are given
    if args.barcode_kit:
        barcode_kit = args.barcode_kit.lower()
        check_guppy_version()
        if args.barcode_kit not in barcode_choices:
            sys.exit('Error: valid --barcodes choices are: {}'.format(join_with_or(barcode_choices)))

    return basecalling_model, barcode_kit

def get_guppy_basecalling_command(input_dir, save_dir, basecalling_model, resume, cpu):
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
    
    print(listToString(basecalling_command))
    return basecalling_command


def get_guppy_barcoder_command(input_dir, save_dir, barcode_kit,resume):
    barcoding_command = ['guppy_barcoder ',
                     '--require_barcodes_both_ends '
                     '--input_path ', input_dir, 
                     '--save_path ', save_dir]
    barcoding_command += BARCODING[barcode_kit]
    if resume:
         barcoding_command += ['--resume']  #add option to resume
    
    print(listToString(barcoding_command))
    return barcoding_command

def get_nextflow_command(demultiplexed_fastq, fast5_pass, sequencing_summary,nf_outdir,run_name,nf_dir_location,conda_location,schemeRepoURL,offline):
    #TODO: add option to specify run folder, cache and medaka options
    logging.info('Running artic guppyplex and artic minion via nextflow pipeline with command: ')
    nextflow_command = ['nextflow run', nf_dir_location,
                     ' -profile conda --nanopolish --cache', conda_location,
                     '--prefix', run_name, 
                     '--basecalled_fastq', demultiplexed_fastq, 
                     '--fast5_pass', fast5_pass,
                     '--sequencing_summary  ', sequencing_summary,
                     '--outdir ', nf_outdir]
    if offline:
        nextflow_command += ['--schemeRepoURL ', schemeRepoURL]  #add option for offline running
    #nextflow_command +=['2>&1 artic_log.txt']

    print(listToString(nextflow_command))
    return nextflow_command

def get_pangolin_command(consensus_file,pangolin_outdir,number_CPUs,offline):
    logging.info('Running pangolin with command: ')
    pangolin_command = ['bash -c "source activate pangolin ; ',]
    if not offline:
        pangolin_command +=[' pangolin --update ; ',] #This takes approx 4 minutes, so consider turning off.
    pangolin_command += [' pangolin ', consensus_file,
                        ' --outdir ', pangolin_outdir,
                        ' --threads ',number_CPUs,
                        ' &>> pangolin_log.txt " ; mv pangolin_log.txt ',pangolin_outdir]

    print(combineCommand(pangolin_command))
    return pangolin_command

def get_nextclade_command(run_name,consensus_dir,nextclade_outdir,offline,dry_run):
    consensus_base=(run_name+'_sequences.fasta')
    if not dry_run or not nextclade_outdir.is_dir():
        os.mkdir(nextclade_outdir) #TODO: ADD TRY
    logging.info('Running nextclade with command: ')
    nextclade_command = []
    if not offline:
        nextclade_command = ['docker pull nextstrain/nextclade ;']  #add option for offline running
    nextclade_command += ['docker run --rm -u 1000' #Note for some systems this is 1000, others 1001
                     ' --volume="',consensus_dir, 
                     ':/seq"  neherlab/nextclade nextclade --input-fasta \'/seq/',consensus_base, 
                     '\' --output-csv \'/seq/',consensus_base,'_nextclade.csv\' ; '
                     'mv ',consensus_dir+consensus_base,'_nextclade.csv ',nextclade_outdir,
                     '&>> ',nextclade_outdir,'nextclade_log.txt ']
    
    print(combineCommand(nextclade_command))
    return nextclade_command

def copy_to_consensus(consensus_dir, artic_outdir, run_name):
    consensus_dir_path=Path(consensus_dir)
    if not consensus_dir_path.is_dir():
        os.mkdir(consensus_dir) #TODO: ADD TRY
    outfilename=(consensus_dir+run_name+'_sequences.fasta')
    with open(outfilename, 'w') as outfile:
        for dirpath, dirs, files in os.walk(artic_outdir):
            for filename in files:
                if filename.endswith('.consensus.fasta'):
                    copyfile(dirpath+'/'+filename, consensus_dir+filename) #TODO: Add sample_name to file?
                    with open(consensus_dir+filename, 'r') as readfile:
                        outfile.write(readfile.read() + "\n\n") 
    consensus_file=outfilename
    return consensus_file

def generate_qc_report(run_name,artic_qc,nextclade_outfile,pangolin_outfile,sample_names,final_report_name, number_CPUs):
    logging.info('Generating run report with QC, lineages and list of mutations')
    #TODO: Add check that files actually exist and that sample_names have correct format:
    print("Found the following files to generate a report for run: " + run_name)
    print("Artic QC file: "+artic_qc)
    print("Nextclade file: "+nextclade_outfile)
    print("Pangolin file: "+pangolin_outfile)
    sample_path = os.path.abspath(sample_names)
    print("Sample names: "+sample_path) #TODO:ignore if header ##TODO: Get full path

    # #Check correct format
    # with open(sample_names, 'r') as f:
    #     print f.readline()

    #Read in files to dataframes
    artic_df = pd.read_csv(artic_qc, sep=',', header=0, encoding='utf8', engine='python')
    nclade_df = pd.read_csv(nextclade_outfile, sep=';', header=0, encoding='utf8', engine='python')
    pangolin_df = pd.read_csv(pangolin_outfile, sep=',', header=0, encoding='utf8', engine='python')
    sample_df = pd.read_csv(sample_names, sep=',', header=0, encoding='utf8', engine='python')

    #TODO: Check that sample_df has header!!
    #If sample_df has "NB[0-9][0-9]" or "Barcode[0-9][0-9]" format, change to "barcode[0-9][0-9]" format
    sample_df['barcode']=sample_df['barcode'].str.replace('NB', 'barcode')
    sample_df['barcode']=sample_df['barcode'].str.replace('Barcode', 'barcode')

    artic_df['run_barcode'] = artic_df.loc[:, 'sample_name']
    nclade_df['run_barcode_artic_nanop'] = nclade_df.loc[:, 'seqName']
    pangolin_df['run_barcode_artic_nanop'] = pangolin_df.loc[:, 'taxon']

    #Get uniform "barcode" and "run" column for each of the dataframes so they can be merged.
    artic_df[['run', 'barcode']] =  artic_df.run_barcode.str.rsplit('_', 1, expand=True).rename(lambda x: f'col{x + 1}', axis=1)
    nclade_df[['run_barcode', 'ARTIC','nanopolish']] = nclade_df['run_barcode_artic_nanop'].str.split('/', 3, expand=True)
    nclade_df[['run', 'barcode']] =  nclade_df.run_barcode.str.rsplit('_', 1, expand=True).rename(lambda x: f'col{x + 1}', axis=1)
    pangolin_df[['run_barcode', 'ARTIC','nanopolish']] = pangolin_df['run_barcode_artic_nanop'].str.split('/', 3, expand=True)
    pangolin_df[['run', 'barcode']] =  pangolin_df.run_barcode.str.rsplit('_', 1, expand=True).rename(lambda x: f'col{x + 1}', axis=1)

    #Merge dataframes
    data_frames = [sample_df, artic_df, pangolin_df, nclade_df]
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['barcode'], how='outer'), data_frames)

    #Create new dataframe with key information
    #From the merged dataframes, take the information you want to be first (rename some of tjen)
    df_main_results = df_merged[['sample_name_x','run_x','barcode','qc_pass','lineage','clade','totalAminoacidSubstitutions']]
    #Extract what you want from the three dataframes and merge as one

    #Merge the new dataframe and the three original dataframes
    final_data_frames = [df_main_results, df_merged]
    df_final_report = reduce(lambda  left,right: pd.merge(left,right,on=['sample_name_x'], how='outer'), final_data_frames)
    #Deselect unneccessary columns
    df_final_report = df_final_report.drop(['run_barcode_artic_nanop_y','run_barcode','ARTIC_y','nanopolish_y','run','run_barcode_artic_nanop_x','run_barcode_y','ARTIC_x','nanopolish_x','run_y','run_barcode_x','run_x_y','barcode_y'], axis=1)
    df_final_report.insert(7, 'artic_QC', 'artic_QC:')
    df_final_report.insert(18, 'pangolin_report', 'pangolin_report:')
    df_final_report.insert(24, 'nextclade_report', 'nextclade_report:')

    df_final_report.rename({'sample_name_x': 'Sample_name', 'run_x_x': 'Run', 'barcode_x': 'Barcode', 'qc_pass_x': 'QC_status', 'lineage_x': 'pangolin_lineage', 'clade_x': 'nextstrain_clade', 'totalAminoacidSubstitutions_x': 'TotalAminoacidSubstitutions', 'sample_name_y': 'sample_name', 'qc_pass_y': 'qc_pass', 'lineage_y': 'lineage', 'clade_y': 'clade', 'totalAminoacidSubstitutions_y': 'totalAminoacidSubstitutions', 'qc.overallStatus': 'qc_overallStatus'}, axis=1, inplace=True)
    
    #Fill any empty rows in the Run col with run_name
    df_final_report.Run = df_final_report.Run.fillna(run_name)

    #Replace TRUE with PASS and FALSE or . with FAIL
    os.environ['NUMEXPR_MAX_THREADS'] = number_CPUs
    df_final_report.QC_status = df_final_report.QC_status.fillna('FAIL')
    mask = df_final_report.applymap(type) != bool
    d = {True: 'PASS', False: 'FAIL'}
    df_final_report = df_final_report.where(mask, df_final_report.replace(d))

    #Add WARN if nextclade report's overall QC is bad:
    df_final_report.loc[df_final_report.qc_overallStatus == 'bad', 'QC_status'] = "WARN"
    df_final_report.loc[df_final_report.qc_pass == 'FAIL', 'QC_status'] = "FAIL"

    #Write to outputfile
    pd.DataFrame.to_csv(df_final_report, final_report_name, sep=',', na_rep='.', index=False)    
    print("Run report written to file: " + final_report_name)

def move_input_files(full_path,raw_data_path,fast5_pass_path,fastq_pass_path,fastq_pass_dem_path):
    #Make 001_rawDAta if not exists
    if not os.path.isdir(raw_data_path):
        os.mkdir(raw_data_path)
    raw_data_path
    #TODO: shutil.move(source,dest) takes too long. Temp solution: Running via bash:
    #Move fast5_pass to 001_rawData
    if not os.path.isdir(os.path.join(raw_data_path, 'fast5_pass')) and os.path.isdir(os.path.join(full_path, 'fast5_pass')):
        source = os.path.join(full_path, 'fast5_pass')
        dest = os.path.join(raw_data_path, 'fast5_pass')
        move_fast5=['mv', source,' ', dest]
        run_command([listToString(move_fast5)], shell=True)

    #Move fastq_pass if it is in input dir
    if not os.path.isdir(os.path.join(raw_data_path, 'fastq_pass')) and os.path.isdir(os.path.join(full_path, 'fastq_pass')):
        source = os.path.join(full_path, 'fastq_pass')
        dest = os.path.join(raw_data_path, 'fastq_pass')
        move_fastq=['mv', source, ' ', dest]
        run_command([listToString(move_fastq)], shell=True)
    #Delete the "work" directory that nextflow leaves behind
    dirpath = os.path.join(full_path, 'work')
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    dirpath = os.path.join(full_path, '.nextflow')
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)

#main
def main():    
    args = parse_args()
    start_time = time.time()
    now = datetime.datetime.now()    
    todays_date = now.strftime('%Y-%m-%d_%H-%M-%S')

    ## Set up log to stdout
    logfile = None
    logging.basicConfig(filename=logfile,level=logging.DEBUG,filemode='w',format='%(asctime)s %(message)s',datefmt='%m-%d-%Y %H:%M:%S')
    logging.info('Running susCovONT v.1.0.0') #Print program version
    logging.info('command line: {0}'.format(' '.join(sys.argv))) #Print input command

    
    #Check with user that input is correct before proceeding. Can be skipped with 'yes |' in beginning of command
    #Check input arguments from config and command and set up environment
    nf_dir_location, conda_location, number_CPUs, schemeRepoURL = set_config_variables()

    full_path, fast5_pass_path, fastq_pass_path, fastq_pass_dem_path, sequencing_summary_path, basecalling_model, barcode_kit, run_name, outdir, nf_outdir, artic_outdir, consensus_dir, pangolin_outdir, nextclade_outdir, nextclade_outfile, pangolin_outfile, artic_qc, sample_name, final_report_name, raw_data_path =check_arguments(args)

    fast5_pass_path, fastq_pass_path, fastq_pass_dem_path, sequencing_summary, print_len_fast5 = check_input(args)


    if not args.dry_run:
        ##Guppy basecalling
        #TODO: add check if basecalling has already been performed
        if args.basecalling_model:
            basecalling_command=(get_guppy_basecalling_command(fast5_pass_path, fastq_pass_path, basecalling_model, args.guppy_resume_basecalling, args.guppy_use_cpu))
            run_command([listToString(basecalling_command)], shell=True)
        
        ##Guppy demultiplexing
        if args.basecalling_model or args.barcode_kit:
            demultiplexing_command=(get_guppy_barcoder_command(fastq_pass_path, fastq_pass_dem_path, barcode_kit,args.guppy_resume_basecalling,))
            run_command([listToString(demultiplexing_command)], shell=True)
        
        #Artic guppyplex and artic minion via PHW's nextflow pipeline
        if not args.generate_report_only and not args.no_artic:
            nextflow_command=(get_nextflow_command(fastq_pass_dem_path, fast5_pass_path, sequencing_summary,nf_outdir,run_name,nf_dir_location,conda_location, schemeRepoURL,args.offline))
            run_command([listToString(nextflow_command)], shell=True)
        consensus_file = copy_to_consensus(consensus_dir,artic_outdir,run_name) ##TODO:FIX mkdir

        ##Pangolin lineage assignment - this worksish  (must add consensus_dir)
        if not args.generate_report_only:
            pangolin_command=(get_pangolin_command(consensus_file,pangolin_outdir,number_CPUs,args.offline))
            run_command([combineCommand(pangolin_command)], shell=True)

        ##Nextclade lineage assignment and substitutions
        if not args.generate_report_only: 
            nextclade_command=(get_nextclade_command(run_name,consensus_dir,nextclade_outdir,args.offline,args.dry_run))
            run_command([combineCommand(nextclade_command)], shell=True)
            
        #Generate run report
        generate_qc_report(run_name,artic_qc,nextclade_outfile,pangolin_outfile,args.sample_names,final_report_name,number_CPUs)
    
        #Move input files to 001_rawData directory
        if not args.no_move_files:
            finish_command=(move_input_files(full_path,raw_data_path,fast5_pass_path,fastq_pass_path,fastq_pass_dem_path))

        #Pipeline complete
        total_time = time.time() - start_time
        time_mins = float(total_time) / 60
        logging.info('susCovONT pipeline finished in ' + str(time_mins) + ' mins.')

    #Dry run
    if args.dry_run:
        ##Guppy basecalling
        if args.basecalling_model:
            basecalling_command=(get_guppy_basecalling_command(fast5_pass_path, fastq_pass_path, basecalling_model, args.guppy_resume_basecalling, args.guppy_use_cpu))
        ##Guppy demultiplexing
        if args.basecalling_model or args.barcode_kit:
            demultiplexing_command=(get_guppy_barcoder_command(fastq_pass_path, fastq_pass_dem_path, barcode_kit,args.guppy_resume_basecalling,))
        #Artic guppyplex and artic minion via PHW's nextflow pipeline
        if not args.generate_report_only and not args.no_artic:
            nextflow_command=(get_nextflow_command(fastq_pass_dem_path, fast5_pass_path, sequencing_summary,nf_outdir,run_name,nf_dir_location,conda_location, schemeRepoURL,args.offline))
        consensus_file = "" #copy_to_consensus(consensus_dir,artic_outdir,run_name)
        ##Pangolin lineage assignment - this worksish  (must add consensus_dir)
        if not args.generate_report_only:
            pangolin_command=(get_pangolin_command(consensus_file,pangolin_outdir,number_CPUs,args.offline))
        ##Nextclade lineage assignment and substitutions
        if not args.generate_report_only: 
            nextclade_command=(get_nextclade_command(run_name,consensus_dir,nextclade_outdir,args.offline,args.dry_run))

    #Pipeline complete
    total_time = time.time() - start_time
    time_mins = float(total_time) / 60
    logging.info('susCovONT pipeline finished in ' + str(time_mins) + ' mins.')

if __name__ == '__main__':
    main()

##TODO:
#TODO: Check if Nextflow has already been run (005_nextclade + "nextflowSucess.txt") and you only want to generate report
#TODO:Check barcode - sample doc
#TODO: Add check for outdir if different than input dir
