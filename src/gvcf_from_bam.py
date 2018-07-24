import subprocess
import glob
import argparse 
import os

def get_args():
    parser = argparse.ArgumentParser(description='given a list of bam files, get the corresponding gvcfs (one by file) and then a poblational gvcf. This script take advantage of gatk docker')
    parser.add_argument('-i','--input_dir',required=True)
    parser.add_argument('-o','--output_dir',required=True) 
    parser.add_argument('-R','--ReferenceFile',required=False, default = '/home/ariel/Projects/BIA/VCFs/data/hgref/human_g1k_v37_decoy.fasta')
    parser.add_argument('-l', '--logfile',required = False, default = './out.log')
    # parse args
    args = parser.parse_args()
    return args

def create_output_dirs(opath):
    if not os.path.exists(opath):
        os.mkdir(opath)

def get_inputfiles(input_dir):
    return glob.glob('%s*.bam'%input_dir)

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



def gatk_docker_get_cmd(basename_file_input, input_dir, output_dir, ReferenceFile, image = 'broadinstitute/gatk'):
    # Reference genonme manipulation and mounting
    ReferenceDir = os.path.dirname(ReferenceFile) # get reference path. It must contain reference.fastafai and reference.dict files. 
    ReferenceFile_basename = os.path.basename(ReferenceFile)


    #Create Docker local paths
    docker_input_data = '/gatk/input_data/'
    docker_dir_output = '/gatk/output/'
    docker_ReferenceDir = '/gatk/refDir/'
    docker_ReferenceFile = '/gatk/refDir/%s'%ReferenceFile_basename 


    ## mounting commands for docker
    mount_input_dir = os.path.abspath(input_dir)+':'+docker_input_data
    mount_output_dir = os.path.abspath(output_dir)+':'+docker_dir_output
    mount_ref_file = os.path.abspath(ReferenceFile)+':'+docker_ReferenceFile
    mount_ref_dir = ReferenceDir+':'+docker_ReferenceDir  

    cmd = ['docker', 'run',
    '-v',mount_input_dir,
    '-v',mount_output_dir,
    '-v',mount_ref_file,
    '-v',mount_ref_dir,
    '-t',image,
    'gatk','--java-options',"-Xmx4g",'HaplotypeCaller',
    '-R',docker_ReferenceFile,
    '-I',docker_input_data + basename_file_input,
    '-O',docker_dir_output + basename_file_input +'.g.vcf',
    '-ERC', 'GVCF'] 

    return cmd



def main():
    args = get_args()
    input_dir, output_dir, ReferenceFile, logfile = args.input_dir, args.output_dir, args.ReferenceFile, args.logfile
    
    create_output_dirs(output_dir) # no es necesario porque el monatje al docker te lo crea si no existe
    inputfiles = get_inputfiles(input_dir)


    for ifile in inputfiles:
        basename_file_input = os.path.basename(ifile)
        
        #create docker call
        cmd = gatk_docker_get_cmd(basename_file_input,input_dir, output_dir, ReferenceFile, image = 'broadinstitute/gatk')
    
        #call docker
        popen_with_logging(cmd,logfile=output_dir + basename_file_input +'_log.out')



if __name__ == '__main__':
    main()