import subprocess
import glob
import argparse 
import os

def get_args():
    parser = argparse.ArgumentParser(description='given a list of bam files, get the corresponding gvcfs (one by file) and then a poblational gvcf. This script take advantage of gatk docker')
    parser.add_argument('-i','--input_files',required=True)
    parser.add_argument('-o','--output_dir',required=True) 
    parser.add_argument('-R','--ReferenceFile',required=False, default = '/home/ariel/Projects/BIA/VCFs/data/hgref/human_g1k_v37_decoy.fasta')
    parser.add_argument('-l', '--logfile',required = False, default = './log.out')
    parser.add_argument('-m', '--mode',required = False, default = 'vcf',help = ' Either "vcf" or "gvcf" ')
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

def get_output_file(input_file, output_dir, mode):
    if mode == 'vcf':
        extension = '.vcf'
    elif mode == 'gvcf':
        extension = '.g.vcf'
    else: 
        print 'mode must be either, vcf or gvcf. EXIT'
        exit
    basename_input_file = os.path.basename(input_file)
    return output_dir + basename_input_file + extension

def gatk_docker_get_cmd(input_file, output_file, ReferenceFile, image = 'broadinstitute/gatk',mode = 'vcf',memory ='4g'):

    # output file extension
    if mode == 'vcf':
        extension = '.vcf'
    elif mode == 'gvcf':
        extension = '.g.vcf'
    else: 
        print 'mode must be either, vcf or gvcf. EXIT'
        exit


    # Reference genonme manipulation and mounting
    ReferenceDir = os.path.dirname(ReferenceFile) # get reference path. It must contain reference.fastafai and reference.dict files. 
    ReferenceFile_basename = os.path.basename(ReferenceFile)
    output_dir = os.path.dirname(output_file)


    #Create Docker local paths
    basename_input_file = os.path.basename(input_file)
    docker_input_file = '/gatk/%s'%basename_input_file
    docker_input_file_index = '/gatk/%s.bai'%basename_input_file

    docker_dir_output = '/gatk/output/'
    docker_ReferenceDir = '/gatk/refDir/'
    docker_ReferenceFile = '/gatk/refDir/%s'%ReferenceFile_basename 


    ## mounting commands for docker
    mount_input_file = os.path.abspath(input_file)+':'+docker_input_file
    mount_output_dir = os.path.abspath(output_dir)+':'+docker_dir_output
    mount_ref_file = os.path.abspath(ReferenceFile)+':'+docker_ReferenceFile
    mount_ref_dir = ReferenceDir+':'+docker_ReferenceDir  

    mount_input_file_index = os.path.abspath(input_file)+'.bai'+':'+docker_input_file_index

    cmd = ['docker', 'run',
    '-v',mount_input_file,
    '-v',mount_input_file_index, 
    '-v',mount_output_dir,
    '-v',mount_ref_file,
    '-v',mount_ref_dir,
    '-t',image,
    'gatk','--java-options',"-Xmx%s"%memory,'HaplotypeCaller',
    '-R',docker_ReferenceFile,
    '-I',docker_input_file,
    '-O',docker_dir_output + basename_input_file + extension]

    if mode == 'gvcf':
        cmd = cmd +['-ERC', 'GVCF'] 

    return cmd



def main():
    args = get_args()
    input_files, output_dir, ReferenceFile, logfile, mode, memory, overwrite = args.input_files, args.output_dir, args.ReferenceFile, \
    args.logfile, args.mode, args.memory, args.overwrite

    create_output_dirs(output_dir) # no es necesario porque el monatje al docker te lo crea si no existe
    inputfiles = glob.glob(input_files)

    print '%s mode activated\n'%mode
    filestorun = "\n"+"\n".join(inputfiles)
    print 'you will run gatk on files: %s \n'%filestorun

    for ifile in inputfiles:
        basename_input_file = os.path.basename(ifile)
        print 'runing %s\n'%basename_input_file
	
        output_file = get_output_file(input_file=ifile,output_dir = output_dir,mode = mode)
        
        if not overwrite:
            if os.path.isfile(output_file):
                print 'file %s already exists, \n switching to NEXT file'%output_file
                continue            

        #create docker call
        cmd = gatk_docker_get_cmd(input_file=ifile, output_file = output_file, ReferenceFile = ReferenceFile, image = 'broadinstitute/gatk',mode = mode,memory=memory)
    	print ' '.join(cmd)
        #call docker
        if logfile is not './log.out':
            logFILE = output_dir + basename_input_file +'.log.out'
        popen_with_logging(cmd,logfile=logFILE)



if __name__ == '__main__':
    main()
