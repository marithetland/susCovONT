#!/usr/bin/env python3
""" marit.hetland@outlook.com
github marithetland, February 2021
This script will basecall and demultiplex
fast5 files from an ONT sequencing run before
running the ARTIC minion pipeline via PHW's Nextflow pipeline. 
Currently tested with versions
Guppy basecalling software: Version 4.4.1
Nextflow v.20.10.0
Artic v1.2.1
ncov nextflow pipeline last updated 23.12.2020 (https://github.com/connor-lab/ncov2019-artic-nf)
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
from functools import reduce
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
from configparser import ConfigParser


#Options for basecalling and barcoding ##Borrowed from Ryan Wick's https://github.com/rrwick/MinION-desktop/blob/master/basecall.py
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


###TODO:s
#Order defs logically (just for fun)
#Add the FASTQ and FAST5 checks in one def , warning if 
###

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
    input_args.add_argument('-s', '--sample_names', type=pathlib.Path, required=True, help='Provide a comma-separated list showing which barcode corresponds to which sample (for final report)')
    
    #Options to perform basecalling and/or demultiplexing before running artic guppyplex and minion
    optional_args.add_argument('--basecall', action='store_true', required=False, help='Perform basecalling before running artic pipeline. Default: off.')
    optional_args.add_argument('--demultiplex', action='store_true', required=False, help='Perform demultiplexing before running artic pipeline. Default: off.')
    optional_args.add_argument('--generate_report_only', action='store_true', required=False, help='Do not run any tools, just (re)generate output report. Default: off.')

    #Only needed if args.basecalling
    optional_args.add_argument('-b', '--basecalling_model', type=str, required=False, choices=["r9.4_fast","r9.4_hac","r10_fast","r10_hac"], help='Indicate which basecalling mode to use. In most cases you want to use a HAC option. This flag is required with --basecalling')
    optional_args.add_argument('--resume_basecalling', action='store_true', required=False, help='This flag can be used with --basecalling to resume an interrupted basecalling run. Default: off.')
    optional_args.add_argument('--cpu', action='store_true', required=False, help='This flag can be used with --basecalling to run on CPU instead of GPU. Will use 4 threads and 6 callers. Default: GPU -auto x.')
    #Only needed if args.demultiplexing
    optional_args.add_argument('-k', '--barcode_kit', type=str, required=False, choices=["none","native_1-12","native_13-24","native_1-24","native_1-96"], help='Indicate which barcode-kits were used, if any. This flag is required with --demultiplexing')

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

def fileCount(path, extension):
    count = 0
    for root, dirs, files in os.walk(path):
        count += sum(f.endswith(extension) for f in files)
    return count

def check_python_version():
    try:
        assert sys.version_info >= (3, 5)
    except AssertionError:
        sys.exit('Error: Python 3.5 or greater is required')

#Only check if guppy_basecalling is being run
def check_guppy_version():
   run_command(['echo -n "guppy_basecaller \t" >> pipeline_versions.txt ; guppy_basecaller --version >> pipeline_versions.txt'], shell=True) ##Check for empty results, skip
   pass

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
    run_command(['echo -n "artic \t" >> pipeline_versions.txt ;~/Programs/conda_for_covid/work/conda/artic-2c6f8ebeb615d37ee3372e543ec21891/bin/artic --version >> pipeline_versions.txt'], shell=True) ##Check for empty results, skip
    pass

#TODO: def check_pangolin_version():
#    run_command(['echo -n "artic \t" >> pipeline_versions.txt ;~/Programs/conda_for_covid/artic-2c6f8ebeb615d37ee3372e543ec21891/bin/artic --version >> pipeline_versions.txt'], shell=True) ##Check for empty results, skip
#    pass


#def check_ncov2019_nextflow_exists():
#    if not os.path.exists(os.path.join(os.getcwd(),fast5_pass)):
#    pass

def check_arguments(args):
    model_choices = list(BASECALLING.keys())
    barcode_choices = list(BARCODING.keys())
    '''First check if basecalling is needed and if correct parameters are given'''
    if args.basecall and (args.basecalling_model is None):
        sys.exit('Error: Flag --basecalling_model must be specified with --basecall. --barcodes choices are: {}'.format(join_with_or(model_choices)))
    if args.basecall and (args.basecalling_model):
        args.basecalling_model = args.basecalling_model.lower()
        check_guppy_version()
        if args.basecalling_model not in model_choices:
            sys.exit('Error: valid --model choices are: {}'.format(join_with_or(model_choices)))

    if args.demultiplex and (args.barcode_kit is None):
        sys.exit('Error: Flag --barcode_kit / -k must be specified with --basecall. --barcodes choices are: {}'.format(join_with_or(barcode_choices)))
    if args.demultiplex and (args.barcode_kit):
        args.barcode_kit = args.barcode_kit.lower()
        check_guppy_version()
        #check_fastq(args)
        if args.barcode_kit not in barcode_choices:
            sys.exit('Error: valid --barcodes choices are: {}'.format(join_with_or(barcode_choices)))

def check_fast5(args): ##This can be done neater
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
                print("Found directory ./001_fast5_pass/ with " + str(len_fast5) + " .fast5 files to analyse")
                raw_fast5s=fast5_pass_alt
            else:
                sys.exit('Error: Found no .fast5 files. Please check that your input directory is correct.')
        else:
            sys.exit('Error: {} is not a directory'.format(fast5_pass))
    if os.path.exists(os.path.join(os.getcwd(),fast5_pass)):
        len_fast5=fileCount(fast5_pass, '.fast5')
        if len_fast5 != 0:
            print("Found directory ./fast5_pass/ with " + str(len_fast5) + " .fast5 files to analyse")
            raw_fast5s=fast5_pass
        else:
            sys.exit('Error: Found no .fast5 files')

    return raw_fast5s

def check_fastq(args): ##This can be done neater
    #Check if the fastq_pass folder is present
    full_path = os.path.abspath(args.input_dir)
    fastq_pass=(full_path+"/fastq_pass/")
    fastq_pass_alt=(full_path+"/002_basecalled/")
    fastq_pass_alt_2=(full_path+"/003_demultiplexed/")
    if os.path.exists(os.path.join(os.getcwd(),fastq_pass)):
        print("Found fastq_pass folder.")
        sequencing_summary=(full_path+"/sequencing_summary*.txt")  #TODO: split run name to the last two "_"s to get this name properly
        basecalled_fastq=fastq_pass
    if not os.path.exists(os.path.join(os.getcwd(),fastq_pass)):
        if os.path.exists(os.path.join(os.getcwd(),fastq_pass_alt)):
            print("Found fastq_pass folder in directory 002_basecalled.")
            sequencing_summary=(fastq_pass_alt+"sequencing_summary.txt")
            if os.path.exists(os.path.join(os.getcwd(),fastq_pass_alt_2)):
                basecalled_fastq=fastq_pass_alt_2
            if not os.path.exists(os.path.join(os.getcwd(),fastq_pass_alt_2)):
                basecalled_fastq=fastq_pass_alt
        else:
            if args.demultiplex:
                print("No fastq_pass folder found. Will perform demultiplexing and store results in ./003_demultiplexed/")
                basecalled_fastq=(full_path+'/002_basecalled/')
                sequencing_summary=(full_path+'/002_basecalled/sequencing_summary.txt')
            if not args.demultiplex:
                sys.exit("No fastq_pass folder found. Make sure the folder 'fastq_pass' or '002_demultiplexed' is in the input_folder, or specifiy --demultiplex to perform demultiplexing of basecalled fast5s." )
    len_fastq=fileCount(basecalled_fastq, '.fastq')
    if len_fastq != 0:
        print("Found directory fastq_pass with " + str(len_fastq) + " .fastq files to analyse")

    return basecalled_fastq, sequencing_summary

def join_with_or(str_list):
    if isinstance(str_list, dict):
        str_list = list(str_list.keys())
    if len(str_list) == 0:
        return ''
    if len(str_list) == 1:
        return str_list[0]
    return ', '.join(str_list[:-1]) + ' or ' + str_list[-1]

    #TODO: Check if Nextflow has already been run (007_nextflow + "nextflowSucess.txt") and you only want to generate report
    #TODO:Check barcode - sample doc
    #TODO: Add check for outdir if different than input dir


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

def get_guppy_barcoder_command(input_dir, save_dir, barcode_kit):
    barcoding_command = ['guppy_barcoder ',
                     '--require_barcodes_both_ends '
                     '--input_path ', input_dir, 
                     '--save_path ', save_dir]
    barcoding_command += BARCODING[barcode_kit]
    #if resume:
    #     guppy_command += ['--resume']  #add option to resume
    
    print(listToString(barcoding_command))
    return barcoding_command

def get_nextflow_command(demultiplexed_fastq, fast5_pass, sequencing_summary,nf_outdir,run_name):
    #TODO: add option to specify run folder, cache and medaka options
    nextflow_command = ['nextflow run ~/Programs/ncov2019-artic-nf -profile conda --nanopolish --cache ~/Programs/conda_for_covid/work/conda ',
                     '--prefix', run_name, 
                     '--basecalled_fastq', demultiplexed_fastq, 
                     '--fast5_pass', fast5_pass,
                     '--sequencing_summary  ', sequencing_summary,
                     '--outdir ', nf_outdir]
    
    print(listToString(nextflow_command))
    return nextflow_command

def get_pangolin_command(consensus_file,pangolin_outdir):
    pangolin_command = ['bash -c "source activate pangolin ; ',
                        #'pangolin --update ',
                        'pangolin ', consensus_file,
                        '--outdir ', pangolin_outdir,
                        '--threads 20 "']
    
    print(listToString(pangolin_command))
    return pangolin_command

def get_nextclade_command(run_name,consensus_dir,nextclade_outdir):
    consensus_base=(run_name+'_sequences.fasta')
    os.mkdir(nextclade_outdir) #TODO: ADD TRY
    nextclade_command = ['docker pull neherlab/nextclade:latest ; ',
                     'docker run -it --rm -u 1001 --volume="',consensus_dir, 
                     ':/seq"  neherlab/nextclade nextclade --input-fasta \'/seq/',consensus_base, 
                     '\' --output-csv \'/seq/',consensus_base,'_nextclade.csv\' ; '
                     'mv ',consensus_dir+consensus_base,'_nextclade.csv ',nextclade_outdir,'']
    
    print(combineCommand(nextclade_command))
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

def copy_to_consensus(consensus_dir, artic_outdir, run_name):
    os.mkdir(consensus_dir) #TODO: ADD TRY
    outfilename=(consensus_dir+run_name+'_sequences.fasta')
    with open(outfilename, 'w') as outfile:
        for dirpath, dirs, files in os.walk(artic_outdir):
            for filename in files:
                if filename.endswith('.consensus.fasta'):
                    os.symlink(dirpath+filename, consensus_dir+filename)
                    with open(consensus_dir+filename, 'r') as readfile:
                        outfile.write(readfile.read() + "\n\n") 
    consensus_file=outfilename
    return consensus_file

def generate_qc_report(run_name,artic_qc,nextclade_outfile,pangolin_outfile,sample_names):
    print("Found the following files to generate a report from: ")
    print("Run name: " +run_name)
    print("Nextclade file: "+nextclade_outfile)
    print("Pangolin file: "+pangolin_outfile)
    print("Artic QC file: "+artic_qc)
    print("Sample names: "+sample_names) #TODO:ignore if header ##TODO: Get full path

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
    df_main_results = df_merged[['sample_name_x','run_x','barcode','qc_pass','lineage','clade','missing','totalMutations','totalAminoacidSubstitutions','totalInsertions','totalAminoacidDeletions','substitutions','aaSubstitutions','insertions','deletions','aaDeletions']]
    #Extract what you want from the three dataframes and merge as one

    #Merge the new dataframe and the three original dataframes
    final_data_frames = [df_main_results, df_merged]
    df_final_report = reduce(lambda  left,right: pd.merge(left,right,on=['sample_name_x'], how='outer'), final_data_frames)
    #Deselect unneccessary columns
    df_final_report = df_final_report.drop(['run_barcode_artic_nanop_y','run_barcode','ARTIC_y','nanopolish_y','run','run_barcode_artic_nanop_x','run_barcode_y','ARTIC_x','nanopolish_x','run_y','run_barcode_x','run_x_y','barcode_y'], axis=1)
    df_final_report.insert(17, 'artic_QC', 'artic_QC:')
    df_final_report.insert(26, 'pangolin_report', 'pangolin_report:')
    df_final_report.insert(34, 'nextclade_report', 'nextclade_report:')

    df_final_report.rename({'sample_name_x': 'Sample_name', 'run_x_x': 'Run', 'barcode_x': 'Barcode', 'qc_pass_x': 'QC_status', 'lineage_x': 'pangolin_lineage', 'clade_x': 'nextstrain_clade', 'missing_x': 'missing_bases', 'totalMutations_x': 'TotalMutations', 'totalAminoacidSubstitutions_x': 'TotalAminoacidSubstitutions', 'totalInsertions_x': 'TotalInsertions', 'totalAminoacidDeletions_x': 'TotalAminoacidDeletions', 'substitutions_x': 'NtSubstitutions', 'aaSubstitutions_x': 'AaSubstitutions','insertions_x': 'NtInsertions', 'deletions_x': 'NtDeletions',  'aaDeletions_x': 'AaDeletions'}, axis=1, inplace=True)

    #Write to outputfile
    final_report_name=(run_name+'_report.txt')
    pd.DataFrame.to_csv(df_final_report, final_report_name, sep=',', na_rep='.', index=False)    
    print("Run report written to file: " + final_report_name)

def get_sample_names(args):
    if not args.sample_names.is_file():
        sys.exit('Error: {} is not a file. Please check your input.'.format(args.sample_names))
    if args.sample_names.is_file():
        sample_names=str(args.sample_names) #TODO: Must check that file exists and has correct header
        print("Your barcodes are specified in: " + sample_names)
    
    #TODO: Check that it is the correct format
    return sample_names

def check_what_to_run(args):

    if args.input_dir:
        full_path = os.path.abspath(args.input_dir)
        run_name = os.path.basename(full_path) 
        print("The input folder is: " + full_path)
        print("The name of your run is: " + run_name)
        outdir=full_path #TODO:Set optional (parental) outdir
        if outdir[-1] != '/':
            outdir = outdir + '/'

        nf_outdir=(outdir+'004_artic_minion/')
        artic_outdir=(outdir+'004_artic_minion/articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish/')
        consensus_dir=(outdir+'005_consensus_fasta/')
        pangolin_outdir=(outdir+'006_pangolin/')
        nextclade_outdir=(outdir+'007_nextclade/')
        nextclade_outfile=(outdir+'007_nextclade/'+run_name+'_sequences.fasta_nextclade.csv')
        pangolin_outfile=(outdir+'006_pangolin/lineage_report.csv')
        artic_qc=(outdir+'004_artic_minion/'+run_name+'.qc.csv')

    if args.basecall and not args.cpu:
        print("Will perform basecalling of fast5 with default GPU")
    if args.basecall and args.cpu:
            print("Will perform basecalling of fast5 in CPU mode (default is GPU)")
    if args.demultiplex:
        print("Will perform demultiplexing of basecalled fast5.")
    if args.generate_report_only:
        print("Will generate final run report from already existing files, no programs will be run.")

    return full_path, run_name, outdir, nf_outdir, artic_outdir, consensus_dir, pangolin_outdir, nextclade_outdir, nextclade_outfile, pangolin_outfile, artic_qc

#main
def main():    
    args = parse_args()
    start_time = time.time()
    now = datetime.datetime.now()    
    todays_date = now.strftime('%Y-%m-%d_%H-%M-%S')

    # ####Sort out config file
    # #Read config.ini file
    # config_object = ConfigParser()
    # config_object.read("config.ini")

    # #Get the password
    # userinfo = config_object["USERINFO"]
    # print("Password is {}".format(userinfo["password"]))
    # ###


    ## Set up log to stdout
    logfile = None
    logging.basicConfig(filename=logfile,level=logging.DEBUG,filemode='w',format='%(asctime)s %(message)s',datefmt='%m-%d-%Y %H:%M:%S')
    logging.info('Running covid-genomics v.1.0.0') #Print program version
    logging.info('command line: {0}'.format(' '.join(sys.argv))) #Print input command

    #raw_fast5s=os.path.abspath(str(args.input_fast5))
    #basecalled_fastq=(outdir+'002_basecalled/') #depends
    #demultiplexed_fastq=(outdir+'003_demultiplexed/') #depends
    #sequencing_summary=(outdir+'002_basecalled/'+'sequencing_summary.txt') #depends on how it was done


    logging.info("##########CHECKPOINT##########")
    full_path, run_name, outdir, nf_outdir, artic_outdir, consensus_dir, pangolin_outdir, nextclade_outdir, nextclade_outfile, pangolin_outfile, artic_qc = check_what_to_run(args)
    raw_fast5=check_fast5(args)
    print(raw_fast5)
    basecalled_fastq,sequencing_summary=check_fastq(args)
    print(basecalled_fastq)
    print(sequencing_summary)
    #Check necessary arguments, set variables and set up output directory
    check_python_version()
    check_guppy_version()
    check_nextflow_version()
    check_artic_version()
    #TODO: # check_ncov2019_nextflow_exists() #check only if not default

    ##Check with user that they are happy with the input for the pipeline
    print("Please check that the specified input below is correct for your run: ")
    #Specify barcode kits and basecaller mode (TODO integrate with def?)
    print("FAST5 files are placed in: ")
    if args.basecall and not args.demultiplex:
        basecaller_mode=args.basecalling_model
        print("Using basecaller mode: " + basecaller_mode)
        print("")
    else:
        print("Basecalled files are placed in: ")
    if args.demultiplex:
        barcode_kit=args.barcode_kit
        print("Specified barcode kit is: " + barcode_kit)
    else:
        print("Basecalled and demultiplexed FASTQ files are placed in: ")

    check_arguments(args)
    print("The input folder is: " + full_path)
    print("The name of your run is: " + run_name)
    sample_names=get_sample_names(args)

    #Check with user that input is correct before proceeding. Can be skipped with 'yes |' in beginning of command
    while True:
        if(yes_or_no('Would you like to run this pipeline? y/n: ')):
            break

    ##Guppy basecalling
    #TODO: add check if basecalling has already been performed
    if args.basecall:
        basecalling_command=(get_guppy_basecalling_command(raw_fast5, basecalled_fastq, basecaller_mode, args.resume_basecalling, args.cpu))
        run_command([listToString(basecalling_command)], shell=True)
    
    ##Guppy demultiplexing
    if args.demultiplex:
        demultiplexing_command=(get_guppy_barcoder_command(basecalled_fastq, demultiplexed_fastq, barcode_kit))
        run_command([listToString(demultiplexing_command)], shell=True)
    
    #Artic guppyplex and artic minion via PHW's nextflow pipeline
    if not args.generate_report_only:
        logging.info('Running artic guppyplex and artic minion')
        nextflow_command=(get_nextflow_command(basecalled_fastq, raw_fast5, sequencing_summary,nf_outdir,run_name))
        run_command([listToString(nextflow_command)], shell=True)
        consensus_file = copy_to_consensus(consensus_dir,artic_outdir,run_name)

    ##Pangolin lineage assignment - this worksish  (must add consensus_dir)
    if not args.generate_report_only:
        logging.info('Running pangolin')
        pangolin_command=(get_pangolin_command(consensus_file,pangolin_outdir))
        run_command([listToString(pangolin_command)], shell=True)

    ##Nextclade lineage assignment and substitutions
    if not args.generate_report_only: 
        logging.info('Running nextclade')
        nextclade_command=(get_nextclade_command(run_name,consensus_dir,nextclade_outdir))
        run_command([combineCommand(nextclade_command)], shell=True)
        
    #Generate run report
    logging.info('Generating run report with QC, lineages and list of mutations')
    qc_command=(generate_qc_report(run_name,artic_qc,nextclade_outfile,pangolin_outfile,sample_names))
   
    #EOF
    total_time = time.time() - start_time
    time_mins = float(total_time) / 60
    logging.info('Covid-genomics pipeline finished in ' + str(time_mins) + ' mins.')

if __name__ == '__main__':
    main()

##TODO:
# Add barkode kits and basecalling model to def and only return/include if it exists.

##How to transfer files to P (and how the file structure should be)
##Easier - how to transfer files from GridION to P
## Note: Add own version of qc.py to run?

##Offline mode?

##Config file:
# Placement of ncov_2019 thingy
# 

#update.py
#Update conda artic
#Update ncov?
#Update pangolin
#Update docker nextclade

##LEgg til sjhekk for NK
