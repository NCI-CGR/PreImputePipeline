#!/usr/bin/env python3

import sys
import math
import os
import subprocess
import shutil
import argparse
import time
import glob

def makeQsub(qsubFile, qsubText):
    '''
    Make a qsub file in the qsubDir
    '''
    with open(qsubFile, 'w') as output:
        output.write('#!/bin/bash\n\n')
        output.write(qsubText)


def runQsub(qsubFile, proj, queue):
    '''
    (str, int) -> None
    '''
    qsubHeader = qsubFile[:-3]
    user = subprocess.check_output('whoami').decode("utf-8").rstrip('\n')
    qsubCall = ['sbatch', '--partition=defq','-o', qsubHeader + '.stdout', '-e', qsubHeader + '.stderr', '--job-name=PreImputeQC.' + proj, qsubFile]
    retcode = subprocess.call(qsubCall)


def getNumSamps(pathToPlink):
    samps = 0
    with open(pathToPlink + '.fam') as f:
        for line in f:
            samps += 1
    return samps

def getNumSnps(pathToPlink):
    snps = 0
    with open(pathToPlink + '.bim') as f:
        for line in f:
            snps += 1
    return snps


def makeConfig(pathToPlink, projDir, reqScripts, strand_chain_file, ref_file, mind, geno, hwe, maf, pop, numSamps):
    '''
    (str, str, str) -> None
    '''
    paths = os.listdir(projDir)
    if 'config.yaml' in paths:
        start = getStartTime(projDir + '/config.yaml')
    else:
        start = time.ctime()
    with open(projDir + '/config.yaml', 'w') as output:
        output.write('plink_genotype_file: ' + str(pathToPlink) + '\n')
        output.write('req_scripts: ' + str(reqScripts) + '\n')
        output.write('strand_chain_file: ' + str(strand_chain_file) + '\n')
        output.write('ref_file: ' + str(ref_file) + '\n')
        output.write('mind: ' + str(mind) + '\n')
        output.write('geno: ' + str(geno) + '\n')
        output.write('hwe: ' + str(hwe) + '\n')
        output.write('maf: ' + str(maf) + '\n')
        output.write('pop: ' + str(pop) + '\n')
        output.write('num_samples: ' + str(numSamps) + '\n')
        output.write('start_time: ' + start + '\n')


def getStartTime(configFile):
    with open(configFile) as f:
        line = f.readline()
        while not line.startswith('start_time'):
            line = f.readline()
        return line.rstrip('\n').split(': ')[1]
    return time.ctime()


def get_args():
    '''
    return the arguments from parser
    '''
    parser = argparse.ArgumentParser()
    requiredArgs = parser.add_argument_group('Required Arguments')
    requiredWithDefaults = parser.add_argument_group('Required arguments with default settings')
    requiredArgs.add_argument('-p', '--path_to_plink_proj', type=str, required=True, help='REQUIRED. Full path to plink binary file project')
    parser.add_argument('-d', '--directory_for_output', type=str, required=True, help='REQUIRED. Full path to the base directory for the pre-imputation pipeline output')
    requiredArgs.add_argument('--req_scripts', type=str, required=True, help='REQUIRED. Full path to required scripts for updating build and HRC_1000G check')
    requiredArgs.add_argument('--strand_chain_file', type=str, required=False, help='REQUIRED. Full path to Strand and liftOver chain files according to genotyping chips and genome build. Files are downloaded from McCarthy Group Tools and liftOver')
    requiredArgs.add_argument('--ref_file', type=str, required=True, help='REQUIRED. Full Path to HRC reference tab file. For example, HRC.r1-1.GRCh37.wgs.mac5.sites.tab')
    requiredWithDefaults.add_argument('--mind', type=float, default= 0.05, help='REQUIRED. Sample call rate filter.  default= 0.95')
    requiredWithDefaults.add_argument('--geno', type=float, default= 0.05, help='REQUIRED. SNP call rate filter.  default= 0.95')
    parser.add_argument('--maf', type=float, default= 0.01, help='OPTIONAL. MAF threshold for SNPs.  default= 0.01')
    parser.add_argument('--hwe', type=float, default= 0.000001, help='OPTIONAL. HWP threshold for SNPs.  default= 0.000001')
    parser.add_argument('--pop', type=str, default= 'ALL', help='OPTIONAL. 1000G population Flag.  default= ALL')
    parser.add_argument('-q', '--queue', type=str, default='defq', help='OPTIONAL. Queue to use to submit jobs.')
    parser.add_argument('-u', '--unlock_snakemake', action='store_true', help='OPTIONAL. If pipeline was killed unexpectedly you may need this flag to rerun')
    args = parser.parse_args()
    return args


def main():
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    args = get_args()
    pathToPlink = args.path_to_plink_proj
    if pathToPlink[0] != '/':
        print('-p argument must be full path to plink project.  Relative paths will not work.')
        sys.exit(1)
    proj = os.path.basename(pathToPlink)
    projDir = args.directory_for_output
    if projDir[0] != '/':
        print('-d argument must be full path to working directory.  Relative paths will not work.')
        sys.exit(1)
    paths = os.listdir(projDir)
#    if 'logs' not in paths:
#        os.mkdir(projDir + '/logs')
    if 'modules' not in paths:
        os.mkdir(projDir + '/modules')
    if not args.strand_chain_file:
        StrandOrChain = None
    else:
        StrChaFile = args.strand_chain_file
        if StrChaFile[-7:] == '.strand':
            StrandOrChain = StrChaFile
        elif StrChaFile[-6:] == '.chain':
            StrandOrChain = StrChaFile
        elif plinkFile[-4:]:
            StrandOrChain = StrChaFile
        else:
            print('Unrecognized Strand or liftOver Chain file format.')
            sys.exit(1)
    #shutil.copyfile('/DCEG/CGF/Bioinformatics/Production/Shilpa/Scripts/PreImputePipeline/Snakefile', projDir + '/Snakefile')
    shutil.copy2(scriptDir + '/Snakefile', projDir + '/Snakefile')
    moduleFiles = glob.glob(scriptDir + '/modules/*')
    for f in moduleFiles:
        shutil.copy2(f, projDir + '/modules')
    numSamps = getNumSamps(pathToPlink)
    numSnps = getNumSnps(pathToPlink)
    makeConfig(pathToPlink, projDir, args.req_scripts, args.strand_chain_file, args.ref_file, args.mind, args.geno, args.hwe, args.maf, args.pop, numSamps)
    qsubTxt = 'cd ' + projDir + '\n'
    qsubTxt += 'module load slurm\n'
    qsubTxt += 'module load python3/3.10.2\n'
    qsubTxt += 'DATE=$(date +%y%m%d)\n'
    qsubTxt += 'mkdir -p logs_${DATE}\n'
    qsubTxt += 'snakemake --cores=2 --unlock --configfile config.yaml\n'
    if args.unlock_snakemake:
        qsubTxt += 'snakemake --unlock\n'
    qsubTxt += 'snakemake -pr --cluster "sbatch --partition=defq --cpus-per-task=10 --output=logs_${DATE}/snakejob_%j.out" --keep-going --rerun-incomplete --jobs 300 --latency-wait 120\n'
    makeQsub(projDir + '/PreImputePipeline.sh', qsubTxt)
    runQsub(projDir + '/PreImputePipeline.sh', proj, args.queue)
    #print('Pre-Imputation Pipeline submitted.  You should receive an email when the pipeline starts and when it completes.')
    print(('Your input file has ' + str(numSamps) + ' samples genotyped at ' + str(numSnps) + ' SNPs'))


if __name__ == "__main__":
    main()
