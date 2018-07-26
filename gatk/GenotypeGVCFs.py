import subprocess
import glob
import argparse 
import os

def get_args():
    parser = argparse.ArgumentParser(description='given a list of bam files, get the corresponding gvcfs (one by file) and then a poblational gvcf. This script take advantage of gatk docker')
    parser.add_argument('-i','--input_files',required=True)
    parser.add_argument('-o','--output_dbpath',required=True) 
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


def GenomicsDBImport_cmd(input_file, output_db_path, image = 'broadinstitute/gatk',docker_db_path = '/gatk/output_db_path/'):

    #Create Docker local paths
    basename_input_file = os.path.basename(input_file)
    docker_input_file = '/gatk/%s'%basename_input_file
    docker_input_file_index = '/gatk/%s.bai'%basename_input_file
    


    ## mounting commands for docker
    mount_input_file = os.path.abspath(input_file)+':'+docker_input_file
    mount_output_dir = os.path.abspath(output_db_path)+':'+docker_db_path
    mount_input_file_index = os.path.abspath(input_file)+'.bai'+':'+docker_input_file_index

    
    cmd = ['docker', 'run',
    '-v',mount_input_file,
    '-v',mount_input_file_index, 
    '-v',mount_output_dir,
    '-t',image,
    'gatk','GenomicsDBImport',
    '-V',docker_input_file,
    '--genomicsdb-workspace-path',
    docker_db_path]

    
    return cmd




def main():
    args = get_args()
    input_files, output_dbpath, logfile, overwrite = args.input_files, args.output_dbpath, args.logfile, args.overwrite

#    create_output_dirs(output_) # no es necesario porque el monatje al docker te lo crea si no existe
    

    #if not overwrite:
    #            if os.path.isdir(output_dbpath):
    #                print 'db folder %s already exists, \n switching to NEXT step (genotypig)'%output_dbpath
    #                continue            


    #create docker call
    cmd = GenomicsDBImport_cmd(input_file=input_files, output_db_path = output_dbpath,image = 'broadinstitute/gatk')
    print ' '.join(cmd)
    #call docker
    basename_input_files = os.path.basename(input_files)
    if logfile is not './log.out':
        logFILE = output_dbpath + basename_input_files +'.log.out'
    popen_with_logging(cmd,logfile=logFILE)



if __name__ == '__main__':
    main()
