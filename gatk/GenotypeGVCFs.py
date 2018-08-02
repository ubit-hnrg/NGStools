# -*- coding: utf-8 -*-
import subprocess
import glob
import argparse 
import os
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description='given a list of bam files, get the corresponding gvcfs (one by file) and then a poblational gvcf. This script take advantage of gatk docker')
    parser.add_argument('-s','--sample_file',required=True)
    parser.add_argument('-d','--sample_path',required=True)
    parser.add_argument('-o','--outpath',required=True) 
    parser.add_argument('-R','--ReferenceFile',required=False, default = '/home/ariel/Projects/BIA/VCFs/data/hgref/human_g1k_v37_decoy.fasta')
    parser.add_argument('-l', '--logfile',required = False, default = './log.out')
    parser.add_argument('-M', '--memory',required = False, default = '4g',help = 'memory to be used in gatk docker')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true')
    parser.add_argument('--no-overwrite', dest='overwrite', action='store_false')
    parser.set_defaults(overwrite=False)
    
    # parse args
    args = parser.parse_args()
    return args

def create_output_dirs(opath):
    if not os.path.exists(opath):
        os.mkdir(opath)


import io
import time
import subprocess
import sys
def popen_with_logging(cmd,logfile = 'out.log'):
    with io.open(logfile, 'wb') as writer, io.open(logfile, 'rb', 1) as reader:
        process = subprocess.Popen(cmd, stdout=writer)
        while process.poll() is None:
            sys.stdout.write(reader.read())
            time.sleep(0.5)
        # Read the remaining
        sys.stdout.write(reader.read())
    return None


def GenomicsDBImport_cmd(sample_file, outpath,sample_path, memory = '4g',image = 'broadinstitute/gatk',dbname = 'mydb',chrom = '12'):

    #Create Docker local paths
    docker_input_dir = '/gatk/inputdata/'
    docker_db_path = '/dbpath/'
    docker_sample_file = '/gatk/samples.txt'


    ## mounting commands for docker
    mount_input_dir = os.path.abspath(sample_path)+':'+docker_input_dir
    mount_output_dir = os.path.abspath(outpath+dbname)+':'+docker_db_path     
    mount_sample_file = sample_file+':'+docker_sample_file 

    
    cmd0 = ['docker', 'run',
    '-v',mount_input_dir,
    '-v',mount_output_dir,
    '-v',mount_sample_file]

    task = ['-t',image,
    'gatk','--java-options','"-Xmx%s -Xms%s"'%(memory,memory),
    'GenomicsDBImport',
    '--genomicsdb-workspace-path',docker_db_path + dbname,
    '--batch-size','50',
    '-L',chrom]
    
    ffiles = ['-V ' + docker_input_dir + line.rstrip('\n') for line in open(sample_file)]    

    cmd = cmd0 + task +ffiles
    
    return cmd

def genotype_gvcf_cmd(ReferenceFile,outpath,dbname,sample_file,memory = '4g',image = 'broadinstitute/gatk'):

    ReferenceDir = os.path.dirname(ReferenceFile)
    ReferenceFile_basename = os.path.basename(ReferenceFile)

    docker_db_path = '/dbpath/'
    docker_outdir = '/outdir/'
    docker_ReferenceDir = '/gatk/refDir/'
    docker_ReferenceFile = '/gatk/refDir/%s'%ReferenceFile_basename

    ## mounting commands for docker
    mount_db_path =  os.path.abspath(outpath+dbname+'/'+dbname)+':'+docker_db_path     
    mount_ref_file = os.path.abspath(ReferenceFile)+':'+docker_ReferenceFile
    mount_ref_dir = ReferenceDir+':'+docker_ReferenceDir
    mount_output_dir =  os.path.abspath(outpath)+ ':' + docker_outdir

    cmd0 = ['docker', 'run',
    '-v',mount_db_path,
    '-v',mount_output_dir,
    '-v',mount_ref_file,
    '-v',mount_ref_dir]

    task = ['-t',image,
        'gatk','--java-options','-Xmx%s'%memory,
        'GenotypeGVCFs',
        '-R', docker_ReferenceFile,
        '-V', 'gendb://'+docker_db_path,
        '-O', docker_outdir + dbname +'_genotypeGVCFs.vcf'
    ]
    cmd = cmd0 + task
    return cmd

def write_and_logging(mje,writer,stdout = True):
    if stdout:
        print mje
    writer.write(mje)

def main():
    args = get_args()
    sample_file, outpath, logfile, sample_path , ReferenceFile,memory = args.sample_file, args.outpath, args.logfile,args.sample_path,args.ReferenceFile,args.memory

    create_output_dirs(outpath) # no es necesario porque el monatje al docker te lo crea si no existe
    CHRMS = [str(i) for i in np.arange(22)+1] +['X','Y']

    #create logging file
    basename_sample_file = os.path.basename(sample_file)
    if logfile is not './log.out':
        logFILE = outpath + basename_sample_file +'.log.out'
    
    #reader1 = io.open(logFILE, 'rb', 1)
    writer1 = io.open(logFILE, 'wb')


    ## dbimport requiere ejecutarse sobre un rango. La manera más practica que veo ahora es x cromosoma, pero tambien puede hacerse via un bed file
    for chrm in CHRMS:
        write_and_logging(mje = '\n doing chr %s \n'%chrm, writer = writer1)

        #create docker call for GenomicDBImport
        dbname = 'dbi'+'_chr'+chrm
        cmd = GenomicsDBImport_cmd(sample_file=sample_file,sample_path = sample_path, dbname=dbname ,outpath = outpath,image = 'broadinstitute/gatk',chrom=chrm,memory=memory)
        write_and_logging(mje = '\n'+' '.join(cmd) + '\n', writer = writer1,stdout=False)        

        #call genomicDBImport via docker
        #popen_with_logging(cmd,logfile=logFILE)
        os.system(" ".join(cmd))   ## WARNING : aca lo force via system, sin logfile, porque popen le pasaba mal el path al docker. No logre entender por que.
        

        # call Genotype via docker
        cmd2 = genotype_gvcf_cmd(ReferenceFile=ReferenceFile,outpath = outpath ,dbname = dbname ,sample_file = sample_file,memory=memory)
        if logfile is not './log.out':
            logFILE = outpath + basename_sample_file +'.log.out'

        write_and_logging(mje = '\n'+' '.join(cmd2) +'\n', writer = writer1,stdout = False)

        popen_with_logging(cmd2,logfile=logFILE)
    
    writer1.close()
if __name__ == '__main__':
    main()
