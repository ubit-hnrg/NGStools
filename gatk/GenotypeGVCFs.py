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
    parser.add_argument('-b','--batchsize',default='10',required = False)
    parser.add_argument('-i','--intervalfile',default=None,required = False)
    parser.add_argument('-D','--dbsnp_path',default='/home/bitgenia/bundle/hgref/',required = False)

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


def GenomicsDBImport_cmd(sample_file, outpath,sample_path, memory = '4g',image = 'broadinstitute/gatk',dbname = 'mydb',interval = None,batchsize = '10'):

    #Create Docker local paths
    docker_input_dir = '/gatk/inputdata/'
    docker_db_path = '/dbpath/'
    docker_sample_file = '/gatk/samples.txt'

    if interval is not None:
        baseint = os.path.basename(interval)
        docker_interval_file = '/gatk/%s'%baseint
        mount_interval_file = interval+':'+docker_interval_file    
    else:
        interval = '1'  # run on chrom 1


    ## mounting commands for docker
    mount_input_dir = os.path.abspath(sample_path)+':'+docker_input_dir
    mount_output_dir = os.path.abspath(outpath+dbname)+':'+docker_db_path     
    mount_sample_file = sample_file+':'+docker_sample_file 

    
    cmd0 = ['docker', 'run',
    '-v',mount_input_dir,
    '-v',mount_output_dir,
    '-v',mount_sample_file,
    '-v',mount_interval_file]

    task = ['-t',image,
    'gatk','--java-options','"-Xmx%s -Xms%s"'%(memory,memory),
    'GenomicsDBImport',
    '--genomicsdb-workspace-path',docker_db_path + dbname,
    '--batch-size',batchsize,
    '-L',docker_interval_file]
    
    ffiles = ['-V ' + docker_input_dir + line.rstrip('\n') for line in open(sample_file)]    

    cmd = cmd0 + task +ffiles
    
    return cmd

def genotype_gvcf_cmd(ReferenceFile,outpath,dbname,sample_file,dbsnp_path = '/home/bitgenia/bundle/hgref/',memory = '4g',image = 'broadinstitute/gatk',interval = None):
    dbsnpfile = 'All_20170710.vcf.gz'

    ReferenceDir = os.path.dirname(ReferenceFile)
    ReferenceFile_basename = os.path.basename(ReferenceFile)

    if interval is not None:
        baseint = os.path.basename(interval)
        docker_interval_file = '/gatk/%s'%baseint
        mount_interval_file = interval+':'+docker_interval_file    
 


    docker_db_path = '/dbpath/'
    docker_outdir = '/outdir/'
    docker_ReferenceDir = '/gatk/refDir/'
    docker_ReferenceFile = '/gatk/refDir/%s'%ReferenceFile_basename
    docker_dbsnp_path = '/gatk/dbsnp_path/'

    ## mounting commands for docker
    mount_db_path =  os.path.abspath(outpath+dbname+'/'+dbname)+':'+docker_db_path     
    mount_ref_file = os.path.abspath(ReferenceFile)+':'+docker_ReferenceFile
    mount_ref_dir = ReferenceDir+':'+docker_ReferenceDir
    mount_output_dir =  os.path.abspath(outpath)+ ':' + docker_outdir
    mount_dbsnp_path =  os.path.abspath(dbsnp_path)+ ':' + docker_dbsnp_path

    cmd0 = ['docker', 'run',
    '-v',mount_db_path,
    '-v',mount_output_dir,
    '-v',mount_ref_file,
    '-v',mount_ref_dir,
    '-v',mount_dbsnp_path,
    '-v',mount_interval_file]

    task = ['-t',image,
        'gatk','--java-options','-Xmx%s'%memory,
        'GenotypeGVCFs',
        '-R', docker_ReferenceFile,
        '-V', 'gendb://'+docker_db_path,
        '-O', docker_outdir + dbname +'_genotypeGVCFs.vcf',
        '-L',docker_interval_file,
        '--dbsnp', docker_dbsnp_path + dbsnpfile,
        '-new-qual'
    ]
    cmd = cmd0 + task
    return cmd

def write_and_logging(mje,writer,stdout = True):
    if stdout:
        print mje
    writer.write(mje)

def main():
    args = get_args()
    sample_file, outpath, logfile, sample_path , ReferenceFile,memory,batchsize,intervalfile, dbsnp_path = args.sample_file, args.outpath, args.logfile,args.sample_path,args.ReferenceFile,args.memory, args.batchsize,args.intervalfile, args.dbsnp_path

    create_output_dirs(outpath) # no es necesario porque el monatje al docker te lo crea si no existe
    ###CHRMS = [str(i) for i in np.arange(22)+1] +['X','Y']

    #create logging file
    basename_sample_file = os.path.basename(sample_file)
    if logfile is not './log.out':
        logFILE = outpath + basename_sample_file +'.log.out'
    
    #reader1 = io.open(logFILE, 'rb', 1)
    writer1 = io.open(logFILE, 'wb')


    ## dbimport requiere ejecutarse sobre un rango. La manera más practica que veo ahora es x cromosoma, pero tambien puede hacerse via un bed file

#    for chrm in CHRMS:
    write_and_logging(mje = '\n doing %s \n'%intervalfile, writer = writer1)

    #create docker call for GenomicDBImport
    dbname = 'dbi'+'_'+ os.path.basename(intervalfile).split('.')[0]
    cmd = GenomicsDBImport_cmd(sample_file=sample_file,sample_path = sample_path, dbname=dbname ,outpath = outpath,image = 'broadinstitute/gatk',interval=intervalfile,memory=memory,batchsize=batchsize)
    write_and_logging(mje = '\n'+' '.join(cmd) + '\n', writer = writer1,stdout=False)        

    #call genomicDBImport via docker
    #popen_with_logging(cmd,logfile=logFILE)
    os.system(" ".join(cmd))   ## WARNING : aca lo force via system, sin logfile, porque popen le pasaba mal el path al docker. No logre entender por que.
    

    # call Genotype via docker
    cmd2 = genotype_gvcf_cmd(ReferenceFile=ReferenceFile,outpath = outpath ,dbname = dbname ,sample_file = sample_file,memory=memory,interval = intervalfile,dbsnp_path = dbsnp_path)
    if logfile is not './log.out':
        logFILE = outpath + basename_sample_file +'.log.out'

    write_and_logging(mje = '\n'+' '.join(cmd2) +'\n', writer = writer1,stdout = False)

    popen_with_logging(cmd2,logfile=logFILE)
    
    writer1.close()
if __name__ == '__main__':
    main()
