#!/usr/bin/python

import json
import os
import argparse

##def comillas_dobles(variable):
##    return '"%s"' %variable


parser = argparse.ArgumentParser(prog='json_upgrade.py',description='generate json input for upgrade_stats', usage='%(prog)s  --samples --bams --bams_index  --trim_front --trim_tail --experiment_lib')
parser.add_argument('-sa','--samples_path',help='samples_path from /data/executionHNRG/experiment_name/inputs')
parser.add_argument('-bam','--bams_path',help='bams_path from /data/resultsHNRG/experiment_name/*.bam')
parser.add_argument('-bai','--bams_index_path',help='bams_index_path from /data/resultsHNRG/experiment_name/*.bai')
parser.add_argument('-tf','--trim_front', help='trimming front bases of fastq for fastp_filtering')
parser.add_argument('-tt','--trim_tail', help='trimming tail bases of fastq for fastp_filtering')
parser.add_argument('-e','--experiment_lib', help='experiment_library')

args = parser.parse_args()
samples = args.samples_path
bams = args.bams_path
bams_index = args.bams_index_path
trim_front = args.trim_front
trim_tail = args.trim_tail
experiment_lib = args.experiment_lib


def upgrade_json(samples,bams,bams_index, trim_front,trim_tail,experiment_lib):
    '''
    Genero un diccionario con los inputs del json original cuyos valores se agregan automaticamente
    '''
    experiment_name = os.path.splitext(os.path.basename(samples))[0]
    output_name = 'inputs_upgrade_stats_'+experiment_name
    data={}
    data['##COMMENT0']='INPUTS automaticos para upgrade_stats. UBIT-HNRG Septiembre 2020'
    data['##COMMENT1']='MANDATORY INPUTS'
    data['##COMMENT2']='INPUT-OUTPUT PATHS'
    data['upgrade_statistics.tabulatedSampleFilePaths']=samples
    data['upgrade_statistics.list_bams']=bams
    data['upgrade_statistics.list_bam_index']=bams_index
    data['upgrade_statistics.experiment_name']=experiment_name+'_upgrade'
    data['upgrade_statistics.path_softlink']= '/data/resultsHNRG/'+experiment_name+'/upgrade_stats_'+trim_front+'-'+trim_tail+'/'
    data['##COMMENT3']='TOOLS & GATK'
    data['upgrade_statistics.toolpath']='/home/hnrg/HNRG-pipeline-V0.1/tools/'
    data['upgrade_statistics.gatk_jar']='gatk-package-4.0.8.1-local.jar'
    data['upgrade_statistics.ngs_toolpath']='/home/hnrg/NGStools'
    data['##COMMENT4']='TRIMMEADO'
    data['upgrade_statistics.trim_front_fastp']=trim_front
    data['upgrade_statistics.trim_tail_fastp']=trim_tail
    data['##COMMENT5']='INTERVALO DE CAPUTURA DEL EXPERIMENTO, ExoNES_ENSEMBL, GENERADOR DE INTERVALO'
    data['upgrade_statistics.padding']='100'
    data['upgrade_statistics.merge_tolerance']='200'
    data['upgrade_statistics.intervalo_captura']=experiment_lib
    data['upgrade_statistics.chromosome_length']='/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/chromosome_lengths_GRCh37_MT.txt'
    data['upgrade_statistics.generic_exon_coords']='/home/hnrg/HNRG-pipeline-V0.1/libraries/intervalos/ensembl_canonicos_GRCh37_0based.tsv'
    out = json.dumps(data,indent = 3)
    return out, output_name, experiment_name
    


json_out, out_name, exp_name = upgrade_json(samples,bams,bams_index, trim_front,trim_tail,experiment_lib)
path_out = '/home/hnrg/executionsHNRG/upgrade_stats/'+exp_name+'/inputs/'

#print(path_out+out_name)
#json_out.to_csv(path_out+out_name+'.json')
with open(path_out+out_name+'.json','w') as outfile:
    outfile.write(json_out)


# }
