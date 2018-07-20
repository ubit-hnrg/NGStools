import os
import glob
import argparse
import fnmatch
import subprocess 
pipeline_path =  '~/repos/TBI/NGS/pipeline/lanes/wesPipeline_lanes_hnrgAddapted.bash' ## hardcoded!! 


def parseargs():
    parser = argparse.ArgumentParser(description='generate lane file paths by sample')
    parser.add_argument('--fastq_folder','-f', dest='fastq_folder',action='store',required = True,help='an integer for the accumulator')
    parser.add_argument('--sep','-s',dest='sep',action= 'store',required = False, default =  '_', help =  'Sample separator in Fastqfile names')
    parser.add_argument('--fastq_extension','-e',dest='fastq_extension',action= 'store',required = False, default =  'fastq', help =  'Fastqfile extension')
    parser.add_argument('--pipeline_output_folder', '-p', dest = 'pipeline_output_folder', action = 'store', required = False, default = './', help = 'output_folder ')
    parser.add_argument('--uncompress', dest='uncompress', action='store_true')
    parser.add_argument('--no-uncompress', dest='uncompress', action='store_false')
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    return args

def get_samplenames(files):
    samplenames = [os.path.basename(f).split(sep)[0] for f in files]
    samplenames = list(set(samplenames)) # drop duplicates
    return samplenames
    
def recursive_glob(treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results

def look_for_sample_paths_and_write_to_disk(samplenames,fastq_extension):
    out = []
    for sample in samplenames:
        r1 = recursive_glob(treeroot= fastq_folder, pattern='%s*R1*%s'%(sample,fastq_extension))
        r2 = recursive_glob(treeroot= fastq_folder, pattern='%s*R2*%s'%(sample,fastq_extension))
        file1 = outputfolder+sample+'_R1.txt'
        with open(file1,'wb') as fr1:
            for line in r1:
                fr1.write(line+'\n')
        fr1.close()        
        file2 = outputfolder+sample+'_R2.txt'
        with open(file2,'wb') as fr2:
            for line in r2:
                fr2.write(line+'\n')
        fr2.close()
        out.append({'sample':sample,'fileR1':file1,'fileR2':file2})
    return out

def uncompress_gz(fastq_extension):
    ff = recursive_glob(treeroot= fastq_folder, pattern='*%s.gz'%fastq_extension)
    for f in ff:
        subprocess.call(["gzip", "-dk",f])
    uncompressed = [f.split('gz')[0] for f in ff]
    return(uncompressed)


def create_pipeline_run(results_folder ,sample,fileR1, fileR2, pipeline_path = pipeline_path , domain = 61, bed_option = 26):
    dom = 'domain=%s'%domain
    bed = 'bed=%s'%bed_option
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)        
    l = ' '.join([pipeline_path,'d',sample,fileR1,fileR2,results_folder,dom,bed])+ '\n'
    return(l)

if __name__ == '__main__':
    args = parseargs()
    fastq_folder, sep, fastq_extension, uncompress , pipeline_output_folder = args.fastq_folder, args.sep, args.fastq_extension, args.uncompress, args.pipeline_output_folder
    files = glob.glob(fastq_folder+'*')
    samplenames = get_samplenames(files)

## uncompress fastq files     
    if uncompress:
        print 'uncompressing files'
        uncompressed = uncompress_gz(fastq_extension)
        print 'files uncompressed'


    outputfolder = '%slanes_samplepath/'%pipeline_output_folder
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)        

    
    list_of_paths = look_for_sample_paths_and_write_to_disk(samplenames,fastq_extension) # esta lista contiene 
    
    ###################################
    ## genero file for running pipeline
    to_run_pipeline = []
    for l in list_of_paths:
        if pipeline_output_folder == './':
            pipeline_output_folder = outputfolder
        to_run_pipeline.append(create_pipeline_run(results_folder = pipeline_output_folder ,sample = l['sample'], fileR1=l['fileR1'], fileR2=l['fileR2']))
    ofile = outputfolder + 'to_run_pipeline.sh'
    
    with open(ofile,'wb') as fo:
        fo.write('#!/bin/bash \n' )
        for line in to_run_pipeline:
            fo.write(line)
    fo.close()
    ###################################

    if uncompress:
        ## gaurdo nombres de archivos descomprimidos para despues de correr el pipeline borrarlos y liberar espacio.
        with open(outputfolder+'uncompressed_files.txt','wb') as f:
            for line in uncompressed:
                f.write(line +'\n')
        f.close()

        
    print 'output is located in the path: \n %s'%outputfolder 
            
        